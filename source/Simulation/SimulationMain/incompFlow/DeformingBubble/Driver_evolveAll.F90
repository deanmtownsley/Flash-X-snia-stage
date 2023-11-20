!!****if* source/Simulation/SimulationMain/incompFlow/DeformingBubble/Driver_evolveAll
!!
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!REORDER(4): face[xyz]Data

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
                             Grid_releaseTileIterator, Grid_getCellCoords, &
                             Grid_restrictAllLevels

   use Grid_iterator, ONLY: Grid_iterator_t

   use Grid_tile, ONLY: Grid_tile_t

   use IO_interface, ONLY: IO_output, IO_outputFinal

   use IncompNS_interface, ONLY: IncompNS_indicators, IncompNS_getGridVar, IncompNS_getScalarProp

   use Multiphase_interface, ONLY: Multiphase_advection, Multiphase_redistance, Multiphase_solve, &
                                   Multiphase_reInitGridVars, Multiphase_indicators, Multiphase_getGridVar

   use Multiphase_data, ONLY: mph_lsIt, mph_extpIt

   use Simulation_data, ONLY: sim_reInitFlow, sim_refineMax

   use RuntimeParameters_interface, ONLY: RuntimeParameters_get

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
   logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
   integer :: iVelVar, iPresVar, iDfunVar, iMfluxVar, &
              iHliqVar, iHgasVar, iTempVar, iDivVar
   integer :: iteration
   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   integer :: ii, jj, kk
   integer, dimension(MDIM) :: lo, hi
   real, allocatable, dimension(:) ::xGrid, yGrid, zGrid
   real :: xi, yi, zi
   logical :: gcell = .true.
   real, pointer, dimension(:, :, :, :) :: facexData, faceyData, facezData
   real, parameter :: pi = acos(-1.0)
   real :: del(MDIM)
   logical :: runUnitTest

   call RuntimeParameters_get('dr_runUnitTest', runUnitTest)

   ! Set interpolation values for guard cell
   call Grid_setInterpValsGcell(.TRUE.)

   ! Get grid variables for incompressible Naiver-Stokes
   call IncompNS_getGridVar("FACE_VELOCITY", iVelVar)

   ! Get grid variables for level set distance function
   call Multiphase_getGridVar("CENTER_LEVELSET", iDfunVar)

   ! Fill GuardCells for level set function
   ! This is done for book-keeping purposes during restart
   call Grid_fillGuardCells(CENTER_FACES, ALLDIR)

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
      dr_simTime = dr_simTime+dr_dt
      dr_simGeneration = 0
      !------------------------------------------------------------

      !------------------------------------------------------------
      ! Set analytical values for fluid flow
      ! sim_reInitFlow is a runtime parameter to force computation
      ! of the velocity field every time step
      !
      ! This is done to test the divergence free interpolation of velocity
      ! If sim_reInitFlow = .FALSE. the divergence should be machine zero
      ! after refinement
      if (sim_reInitFlow) then

         ! Loop over iterator
         call Grid_getTileIterator(itor, nodetype=LEAF)
         do while (itor%isValid())
            call itor%currentTile(tileDesc)

            ! Logic to compute analytical solution of velocity
            !---------------------------------------------------------
            nullify (facexData, faceyData, facezData)

            lo = tileDesc%blkLimitsGC(LOW, :)
            hi = tileDesc%blkLimitsGC(HIGH, :)

            call tileDesc%deltas(del)

            allocate (xGrid(lo(IAXIS):hi(IAXIS)+1))
            allocate (yGrid(lo(JAXIS):hi(JAXIS)))

            xGrid = 0.0
            yGrid = 0.0

            call Grid_getCellCoords(IAXIS, FACES, tileDesc%level, lo, hi, xGrid)
            call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, lo, hi, yGrid)

            call tileDesc%getDataPtr(facexData, FACEX)
            do kk = lo(KAXIS), hi(KAXIS)
               do jj = lo(JAXIS), hi(JAXIS)
                  do ii = lo(IAXIS), hi(IAXIS)+1
                     xi = xGrid(ii)
                     yi = yGrid(jj)

                     facexData(VELC_FACE_VAR, ii, jj, kk) = -((sin(pi*xi))**2)* &
                                                            (cos(2*pi*(yi+del(JAXIS)/2))- &
                                                             cos(2*pi*(yi-del(JAXIS)/2)))/(2*pi*del(JAXIS))

                  end do
               end do
            end do
            call tileDesc%releaseDataPtr(facexData, FACEX)
            deallocate (xGrid, yGrid)

            allocate (xGrid(lo(IAXIS):hi(IAXIS)))
            allocate (yGrid(lo(JAXIS):hi(JAXIS)+1))

            xGrid = 0.0
            yGrid = 0.0

            call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xGrid)
            call Grid_getCellCoords(JAXIS, FACES, tileDesc%level, lo, hi, yGrid)

            call tileDesc%getDataPtr(faceyData, FACEY)
            do kk = lo(KAXIS), hi(KAXIS)
               do jj = lo(JAXIS), hi(JAXIS)+1
                  do ii = lo(IAXIS), hi(IAXIS)
                     xi = xGrid(ii)
                     yi = yGrid(jj)

                     faceyData(VELC_FACE_VAR, ii, jj, kk) = ((sin(pi*yi))**2)* &
                                                            (cos(2*pi*(xi+del(IAXIS)/2))- &
                                                             cos(2*pi*(xi-del(IAXIS)/2)))/(2*pi*del(IAXIS))

                  end do
               end do
            end do
            call tileDesc%releaseDataPtr(faceyData, FACEY)
            deallocate (xGrid, yGrid)

#if NDIM == MDIM
            call tileDesc%getDataPtr(facezData, FACEZ)
            facezData(VELC_FACE_VAR, :, :, :) = 0.
            call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif
            !---------------------------------------------------------

            call itor%next()
         end do
         call Grid_releaseTileIterator(itor)
      end if
      !------------------------------------------------------------

      ! Compute divergences for velocities
      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call IncompNS_divergence(tileDesc)
         !---------------------------------------------------------
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
      !------------------------------------------------------------

      ! Call indicators method to show information
      !------------------------------------------------------------
      call IncompNS_indicators()
      !------------------------------------------------------------

      ! Call methods to reset specific grid variables
      ! for Multiphase unit
      !------------------------------------------------------------
      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call Multiphase_reInitGridVars(tileDesc)
         !---------------------------------------------------------
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
      !------------------------------------------------------------

      ! Fill GuardCells for level set function
      gcMask(:) = .FALSE.
      gcMask(iDfunVar) = .TRUE.
      call Grid_fillGuardCells(CENTER, ALLDIR, &
                               maskSize=NUNK_VARS, mask=gcMask)

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

      ! Apply redistancing procedure
      !------------------------------------------------------------
      do iteration = 1, mph_lsIt

         ! Fill GuardCells for level set function
         gcMask(:) = .FALSE.
         gcMask(iDfunVar) = .TRUE.
         call Grid_fillGuardCells(CENTER, ALLDIR, &
                                  maskSize=NUNK_VARS, mask=gcMask)

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

      end do
      !------------------------------------------------------------

      ! Call indicators to show information
      !------------------------------------------------------------
      call Multiphase_indicators()
      !------------------------------------------------------------

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
                     dr_dt, dr_nstep+1, dr_nbegin, endRunPl, PLOTFILE_AND_PARTICLEFILE)
      call Timers_stop("IO_output")

      ! Update grid and notify changes to other units
      call Grid_updateRefinement(dr_nstep, dr_simTime, gridChanged)

      if (gridChanged) then
         dr_simGeneration = dr_simGeneration+1
         if (runUnitTest) call Grid_fillGuardCells(CENTER_FACES, ALLDIR)
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
      call IO_output(dr_simTime, dr_dt, dr_nstep+1, dr_nbegin, endRun, &
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
   call Timers_getSummary(max(0, dr_nstep-dr_nbegin+1))
   call Logfile_stamp("FLASH run complete.", "LOGFILE_END")
   call Logfile_close()

   !-------------------------------------------------------------------------------
   !--uniTest procedures------
   call IncompNS_getScalarProp("Min_Divergence", mindiv)
   call IncompNS_getScalarProp("Max_Divergence", maxdiv)

   if (runUnitTest) then
      temp = dr_globalMe
      do i = 1, 4
         prNum(i) = mod(temp, 10)
         temp = temp/10
      end do
      filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))// &
                 char(48+prNum(2))//char(48+prNum(1))
      open (fileUnit, file=fileName)
      write (fileUnit, '("P",I0)') dr_globalMe
      if (abs(mindiv) .le. 1e-11 .and. abs(maxdiv) .le. 1e-11) then
         write (fileUnit, '(A)') 'SUCCESS all results conformed with expected values.'
      else
         write (fileUnit, '(A)') 'FAILURE'
      end if
      close (fileUnit)
   end if
   !-------------------------------------------------------------------------------

   return

end subroutine Driver_evolveAll
