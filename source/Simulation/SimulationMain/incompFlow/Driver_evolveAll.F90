!!****if* source/Simulation/SimulationMain/incompFlow/Driver_evolveAll
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!! NAME
!!
!!  Driver_evolveAll
!!
!! SYNOPSIS
!!
!!  call Driver_evolveAll()
!!
!! DESCRIPTION
!!
!!  Driver_evolveAll for incompFlow Simulations
!!
!!  DOC: Driver_evolveAll needs more explanation
!!
!! NOTES
!!
!!  variables that begin with "dr_" like, dr_globalMe or dr_dt, dr_beginStep
!!  are stored in the data fortran module for the Driver unit, Driver_data.
!!  The "dr_" is meant to indicate that the variable belongs to the Driver Unit.
!!  all other normally named variables i, j, etc are local variables.
!!
!!
!!***

#include "constants.h"
#include "Simulation.h"

#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_evolveAll()

   use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs, dr_nbegin, &
                          dr_nend, dr_dt, &
                          dr_tmax, dr_simTime, dr_redshift, &
                          dr_nstep, dr_dtOld, dr_dtNew, &
                          dr_simGeneration, &
                          dr_restart
   use Driver_interface, ONLY: Driver_sourceTerms, Driver_computeDt
   use Logfile_interface, ONLY: Logfile_stamp, Logfile_close
   use Timers_interface, ONLY: Timers_start, Timers_stop, &
                               Timers_getSummary
   use Particles_interface, ONLY: Particles_advance, Particles_dump

   use Grid_interface, ONLY: Grid_updateRefinement, Grid_setInterpValsGcell, &
                             Grid_fillGuardCells, Grid_getTileIterator, &
                             Grid_releaseTileIterator, Grid_solvePoisson

   use Grid_iterator, ONLY: Grid_iterator_t

   use Grid_tile, ONLY: Grid_tile_t

   use IO_interface, ONLY: IO_output, IO_outputFinal

   use IncompNS_interface, ONLY: IncompNS_predictor, IncompNS_divergence, &
                                 IncompNS_setupPoisson, IncompNS_corrector, &
                                 IncompNS_indicators, IncompNS_reInitGridVars, &
                                 IncompNS_advection, IncompNS_diffusion, IncompNS_getScalarProp, &
                                 IncompNS_getGridVar

   use Multiphase_interface, ONLY: Multiphase_setFluidProps, Multiphase_setPressureJumps, &
                                   Multiphase_advection, Multiphase_redistance, Multiphase_solve, &
                                   Multiphase_setThermalProps, Multiphase_thermalForcing, Multiphase_velForcing, &
                                   Multiphase_divergence, Multiphase_extrapFluxes, Multiphase_reInitGridVars, &
                                   Multiphase_indicators, Multiphase_setMassFlux, Multiphase_getGridVar

   use HeatAD_interface, ONLY: HeatAD_diffusion, HeatAD_advection, HeatAD_solve, HeatAD_reInitGridVars, &
                               HeatAD_indicators, HeatAD_getGridVar

   use Simulation_interface, ONLY: Simulation_adjustEvolution

   use IncompNS_data, ONLY: ins_predcorrflg, ins_pressureBC_types, ins_pressureBC_values, &
                            ins_poisfact

#ifdef MULTIPHASE_MAIN
   use Multiphase_data, ONLY: mph_lsIt, mph_extpIt, mph_iJumpVar
#endif

   implicit none

   ! for logfile output
   character(len=MAX_STRING_LENGTH), dimension(3, 2) :: strBuff
   character(len=15) :: numToStr

   logical :: gridChanged
   logical :: endRunPl !Should we end our run on this iteration, based on conditions detected by the IO unit?
   logical :: endRun !Should we end our run on this iteration, based on conditions detected by the IO unit?
   logical :: endRunWallClock !Should we end our run on this iteration, based on wall clock time?

   !---------Variables for unitTest-----------
   character(len=20) :: fileName
   integer, parameter        :: fileUnit = 2
   integer, dimension(4) :: prNum
   integer :: temp, i
   real :: mindiv, maxdiv
   logical :: gcMask(NUNK_VARS + NDIM*NFACE_VARS)
   integer :: iVelVar, iPresVar, iDfunVar, iMfluxVar, &
              iHliqVar, iHgasVar, iTempVar, iDivVar
   integer :: iteration
   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   ! Get grid variables for incompressible Naiver-Stokes
   call IncompNS_getGridVar("FACE_VELOCITY", iVelVar)
   call IncompNS_getGridVar("CENTER_PRESSURE", iPresVar)
   call IncompNS_getGridVar("CENTER_DIVERGENCE", iDivVar)

#ifdef HEATAD_MAIN
   ! Get grid variables for heat advection diffusion if
   ! if HeatAD unit is available
   call HeatAD_getGridVar("CENTER_TEMPERATURE", iTempVar)
#endif

#ifdef MULTIPHASE_MAIN
   ! Get grid variables for level set distance function
   ! if Multiphase unit is available
   call Multiphase_getGridVar("CENTER_LEVELSET", iDfunVar)

   ! Fill GuardCells for level set function
   ! This is done for book-keeping purposes during
   ! restart
   gcMask = .FALSE.
   gcMask(iDfunVar) = .TRUE.
   call Grid_fillGuardCells(CENTER, ALLDIR, &
                            maskSize=NUNK_VARS, mask=gcMask)
#endif

#ifdef MULTIPHASE_EVAPORATION
   ! Get grid variables for heat flux and mass flux if
   ! Multiphase evaporation is present
   call Multiphase_getGridVar("CENTER_MASSFLUX", iMfluxVar)
   call Multiphase_getGridVar("CENTER_HFLUX_LIQUID", iHliqVar)
   call Multiphase_getGridVar("CENTER_HFLUX_GAS", iHgasVar)
#endif

   endRunPl = .false.
   endRun = .false.

   call Logfile_stamp('Entering evolution loop', '[Driver_evolveAll]')
   call Timers_start("evolution")

   ! Initial Timestep:
   ! backup needed old
   dr_dtOld = dr_dt

   ! calculate new
   call Driver_computeDt(dr_nbegin, dr_nstep, &
                         dr_simTime, dr_dtOld, dr_dtNew)
   ! store new
   dr_dt = dr_dtNew

   do dr_nstep = dr_nBegin, dr_nend

      if (dr_globalMe == MASTER_PE) then

         write (numToStr(1:), '(I10)') dr_nstep
         write (strBuff(1, 1), "(A)") "n"
         write (strBuff(1, 2), "(A)") trim(adjustl(numToStr))

         write (numToStr(1:), "(1PE12.6)") dr_simTime
         write (strBuff(2, 1), "(A)") "t"
         write (strBuff(2, 2), "(A)") trim(adjustl(numToStr))

         write (numToStr(1:), "(1PE12.6)") dr_dt
         write (strBuff(3, 1), "(A)") "dt"
         write (strBuff(3, 2), "(A)") trim(adjustl(NumToStr))

         call Logfile_stamp(strBuff, 3, 2, "step")
      end if

      !------------------------------------------------------------
      !- Start Physics Sequence
      !------------------------------------------------------------
      dr_simTime = dr_simTime + dr_dt
      dr_simGeneration = 0
      !------------------------------------------------------------

      ! Call methods to reset specific grid variables
      ! at the start of every time-setp
      !------------------------------------------------------------
      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call Multiphase_reInitGridVars(tileDesc)
         call IncompNS_reInitGridVars(tileDesc)
         call HeatAD_reInitGridVars(tileDesc)
         !---------------------------------------------------------
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
      !------------------------------------------------------------

      ! Apply simulation specific modifications
      !------------------------------------------------------------
      call Simulation_adjustEvolution(dr_nstep, dr_dt, dr_simTime)
      !------------------------------------------------------------

#ifdef MULTIPHASE_MAIN
      ! Multiphase advection procedure
      ! Loop over blocks (tiles) and call Multiphase
      ! routines
      !------------------------------------------------------------
      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call Multiphase_advection(tileDesc)
         call Multiphase_solve(tileDesc, dr_dt)
         !---------------------------------------------------------
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
      !------------------------------------------------------------

      ! Fill GuardCells for level set function
      gcMask = .FALSE.
      gcMask(iDfunVar) = .TRUE.
      call Grid_fillGuardCells(CENTER, ALLDIR, &
                               maskSize=NUNK_VARS, mask=gcMask)

      ! Apply redistancing procedure
      !------------------------------------------------------------
      do iteration = 1, mph_lsIt

         ! Loop over blocks (tiles) and call Multiphase
         ! routines
         call Grid_getTileIterator(itor, nodetype=LEAF)
         do while (itor%isValid())
            call itor%currentTile(tileDesc)
            !------------------------------------------------------
            call Multiphase_redistance(tileDesc, iteration)
            !------------------------------------------------------
            call itor%next()
         end do
         call Grid_releaseTileIterator(itor)

         ! Fill GuardCells for level set function
         gcMask = .FALSE.
         gcMask(iDfunVar) = .TRUE.
         call Grid_fillGuardCells(CENTER, ALLDIR, &
                                  maskSize=NUNK_VARS, mask=gcMask)
      end do
      !------------------------------------------------------------

      ! Call indicators to show information
      !------------------------------------------------------------
      call Multiphase_indicators()
      !------------------------------------------------------------

      ! Update fluid properties and pressure jumps
      ! Loop over blocks (tiles)
      !------------------------------------------------------------
      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call Multiphase_setFluidProps(tileDesc)
         call Multiphase_setPressureJumps(tileDesc)
         !---------------------------------------------------------
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
      !------------------------------------------------------------

      ! Fill GuardCells for Pressure Jump
      gcMask = .FALSE.
      gcMask(NUNK_VARS + mph_iJumpVar) = .TRUE.
      gcMask(NUNK_VARS + 1*NFACE_VARS + mph_iJumpVar) = .TRUE.
#if NDIM == 3
      gcMask(NUNK_VARS + 2*NFACE_VARS + mph_iJumpVar) = .TRUE.
#endif
      call Grid_fillGuardCells(CENTER_FACES, ALLDIR, &
                               maskSize=NUNK_VARS + NDIM*NFACE_VARS, mask=gcMask)

#ifdef MULTIPHASE_EVAPORATION
      ! Update thermal properties and apply thermal
      ! forcing for evaporation
      !------------------------------------------------------------
      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call Multiphase_setThermalProps(tileDesc)
         call Multiphase_thermalForcing(tileDesc)
         !---------------------------------------------------------
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
      !------------------------------------------------------------

      ! Perform extrapolation iterations for
      ! heat flux
      !------------------------------------------------------------
      do iteration = 1, mph_extpIt

         call Grid_getTileIterator(itor, nodetype=LEAF)
         do while (itor%isValid())
            call itor%currentTile(tileDesc)
            !------------------------------------------------------
            call Multiphase_extrapFluxes(tileDesc, iteration)
            !------------------------------------------------------
            call itor%next()
         end do
         call Grid_releaseTileIterator(itor)

         ! Fill GuardCells for heat fluxes
         gcMask = .FALSE.
         gcMask(iHliqVar) = .TRUE.
         gcMask(iHGasVar) = .TRUE.
         call Grid_fillGuardCells(CENTER, ALLDIR, &
                                  maskSize=NUNK_VARS, mask=gcMask)
      end do
      !------------------------------------------------------------

      ! Set mass flux from extrapolated heat fluxes
      !------------------------------------------------------------
      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call Multiphase_setMassFlux(tileDesc)
         !---------------------------------------------------------
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
      !------------------------------------------------------------

      ! Fill GuardCells for mass flux
      gcMask = .FALSE.
      gcMask(iMfluxVar) = .TRUE.
      call Grid_fillGuardCells(CENTER, ALLDIR, &
                               maskSize=NUNK_VARS, mask=gcMask)
#endif

#endif

      !------------------------------------------------------------
      call Grid_setInterpValsGcell(.true.)
      !------------------------------------------------------------

      ! Fill GuardCells for velocity
      gcMask = .FALSE.
      gcMask(NUNK_VARS + iVelVar) = .TRUE.
      gcMask(NUNK_VARS + 1*NFACE_VARS + iVelVar) = .TRUE.
#if NDIM == 3
      gcMask(NUNK_VARS + 2*NFACE_VARS + iVelVar) = .TRUE.
#endif
      call Grid_fillGuardCells(CENTER_FACES, ALLDIR, &
                               maskSize=NUNK_VARS + NDIM*NFACE_VARS, mask=gcMask)

      ! Start of fractional-step velocity procedure
      ! Calculate predicted velocity and apply
      ! necessary forcing
      !------------------------------------------------------------
      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call IncompNS_advection(tileDesc)
         call IncompNS_diffusion(tileDesc)
         call IncompNS_predictor(tileDesc, dr_dt)
         call Multiphase_velForcing(tileDesc, dr_dt)
         !---------------------------------------------------------
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
      !------------------------------------------------------------

      ! Fill GuardCells for velocity
      gcMask = .FALSE.
      gcMask(NUNK_VARS + iVelVar) = .TRUE.
      gcMask(NUNK_VARS + 1*NFACE_VARS + iVelVar) = .TRUE.
#if NDIM == 3
      gcMask(NUNK_VARS + 2*NFACE_VARS + iVelVar) = .TRUE.
#endif
      ins_predcorrflg = .true.
      call Grid_fillGuardCells(CENTER_FACES, ALLDIR, &
                               maskSize=NUNK_VARS + NDIM*NFACE_VARS, mask=gcMask)
      ins_predcorrflg = .false.

      ! Calculate divergence of predicted velocity
      !------------------------------------------------------------
      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call IncompNS_divergence(tileDesc)
         call Multiphase_divergence(tileDesc)
         call IncompNS_setupPoisson(tileDesc, dr_dt)
         !---------------------------------------------------------
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
      !------------------------------------------------------------

      ! Solve pressure Poisson equation
      !------------------------------------------------------------
      call Grid_solvePoisson(iSoln=iPresVar, iSrc=iDivVar, &
                             bcTypes=ins_pressureBC_types, &
                             bcValues=ins_pressureBC_values, &
                             poisfact=ins_poisfact)
      !------------------------------------------------------------

      ! Fill GuardCells for pressure
      gcMask = .FALSE.
      gcMask(iPresVar) = .TRUE.
      call Grid_fillGuardCells(CENTER, ALLDIR, &
                               maskSize=NUNK_VARS, mask=gcMask, &
                               selectBlockType=ACTIVE_BLKS)

      ! Final step of fractional step velocity
      ! formulation - calculate corrected velocity
      ! and updated divergence (this should be machine-zero)
      !------------------------------------------------------------
      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call IncompNS_corrector(tileDesc, dr_dt)
         call IncompNS_divergence(tileDesc)
         call Multiphase_divergence(tileDesc)
         !---------------------------------------------------------
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
      !------------------------------------------------------------

      ! Call indicators method to show information
      !------------------------------------------------------------
      call IncompNS_indicators()
      !------------------------------------------------------------

      !------------------------------------------------------------
      call Grid_setInterpValsGcell(.false.)
      !------------------------------------------------------------

#ifdef HEATAD_MAIN
      ! Fill GuardCells for temperature
      gcMask = .FALSE.
      gcMask(iTempVar) = .TRUE.
      call Grid_fillGuardCells(CENTER, ALLDIR, &
                               maskSize=NUNK_VARS, mask=gcMask)

      ! Heat advection diffusion procedure
      !------------------------------------------------------------
      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call HeatAD_diffusion(tileDesc)
         call HeatAD_advection(tileDesc)
         call HeatAD_solve(tileDesc, dr_dt)
         !---------------------------------------------------------
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
      !------------------------------------------------------------

      ! Call indicators
      !------------------------------------------------------------
      call HeatAD_indicators()
      !------------------------------------------------------------
#endif

      !------------------------------------------------------------
      !- End Physics Sequence
      !------------------------------------------------------------
      ! Velocities and Omg to Center variables
      ! In your Simulation Config set REQUIRES physics/IncompNS/IncompNSExtras
      ! Note this will add velocity and vorticity variables to your CENTER data structure.
      ! Average Velocities and Vorticity to cell-centers
      call IncompNS_velomgToCenter()

      !output a plotfile before the grid changes
      call Timers_start("IO_output")

      call IO_output(dr_simTime, &
                     dr_dt, dr_nstep + 1, dr_nbegin, endRunPl, PLOTFILE_AND_PARTICLEFILE)
      call Timers_stop("IO_output")

      ! Update grid and notify changes to other units
      call Grid_updateRefinement(dr_nstep, dr_simTime, gridChanged)

      ! Perform housekeeping after gridChanges
      ! Fill guard cells for new grid
      if (gridChanged) then
         dr_simGeneration = dr_simGeneration + 1

         gcMask = .FALSE.
         gcMask(iDfunVar) = .TRUE.
         call Grid_fillGuardCells(CENTER, ALLDIR, &
                                  maskSize=NUNK_VARS, mask=gcMask)
      end if

      if (dr_globalMe .eq. MASTER_PE) then
         write (*, *) ' '
         write (*, '(I6,A,g16.8,A,g16.8)') dr_nstep, &
            ', TimeStep= ', dr_dt, ', SimTime= ', dr_simTime
      end if

      if (dr_globalMe .eq. MASTER_PE) &
         write (*, *) '###############################################################################'

      ! Compute next step dt:
      ! backup needed old
      dr_dtOld = dr_dt

      ! calculate new
      call Driver_computeDt(dr_nbegin, dr_nstep, &
                            dr_simTime, dr_dtOld, dr_dtNew)
      ! store new
      dr_dt = dr_dtNew

      call Timers_start("io")
      call IO_output(dr_simTime, dr_dt, dr_nstep + 1, dr_nbegin, endRun, &
                     CHECKPOINT_FILE_ONLY)
      call Timers_stop("io")
      endRun = (endRunPl .OR. endRun)

     !!*****************************************************************************
     !!  Evolution Loop -- check termination conditions
     !!*****************************************************************************

      !Exit if a .dump_restart or .kill was found during the last step
      if (endRun) exit

     !! the simulation ends before nend iterations if
     !!  (i)   the simulation time is greater than the maximum time (tmax)
     !!  (ii)  the redshift falls below the minimum redshift
     !!        (also called zfinal)
     !!  (iii) the wall clock time is greater than the maximum
     !!        (wall_clock_time_max)

      if (dr_simTime >= dr_tmax) then
         if (dr_globalMe == MASTER_PE) then
            print *, "exiting: reached max SimTime"
         end if
         exit
      end if

      call dr_wallClockLimitExceeded(endRunWallClock)
      if (endRunWallClock) then
         if (dr_globalMe == MASTER_PE) then
            print *, "exiting: reached max wall clock time"
         end if
         exit
      end if

   end do
   dr_nstep = min(dr_nstep, dr_nend)

  !!******************************************************************************
  !! End of Evolution Loop
  !!******************************************************************************

   call Timers_stop("evolution")
   call Logfile_stamp('Exiting evolution loop', '[Driver_evolveAll]')
   if (.NOT. endRun) call IO_outputFinal()
   call Timers_getSummary(max(0, dr_nstep - dr_nbegin + 1))
   call Logfile_stamp("FLASH run complete.", "LOGFILE_END")
   call Logfile_close()

   !-------------------------------------------------------------------------------
   !--uniTest procedures------
   call IncompNS_getScalarProp("Min_Divergence", mindiv)
   call IncompNS_getScalarProp("Max_Divergence", maxdiv)

   temp = dr_globalMe
   do i = 1, 4
      prNum(i) = mod(temp, 10)
      temp = temp/10
   end do
   filename = "unitTest_"//char(48 + prNum(4))//char(48 + prNum(3))// &
              char(48 + prNum(2))//char(48 + prNum(1))
   open (fileUnit, file=fileName)
   write (fileUnit, '("P",I0)') dr_globalMe
   if (abs(mindiv) .le. 1e-11 .and. abs(maxdiv) .le. 1e-11) then
      write (fileUnit, '(A)') 'SUCCESS all results conformed with expected values.'
   else
      write (fileUnit, '(A)') 'FAILURE'
   end if
   close (fileUnit)
   !-------------------------------------------------------------------------------

   return

end subroutine Driver_evolveAll
