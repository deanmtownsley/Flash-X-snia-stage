!!****if* source/Simulation/SimulationMain/heatDiffusion
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
                             Grid_fillGuardCells

   use IO_interface, ONLY: IO_output, IO_outputFinal

   use HeatAD_interface, ONLY: HeatAD_diffusion, HeatAD_solve, HeatAD_reInitGridVars, &
                               HeatAD_indicators, HeatAD_getGridVar

   use Simulation_interface, ONLY: Simulation_adjustEvolution

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
   integer :: iVelVar, iPresVar, iDfunVar, iMfluxVar, iHliqVar, iHgasVar, iTempVar
   integer :: iteration

   ! Get grid variables for heat diffusion
   call HeatAD_getGridVar("CENTER_TEMPERATURE", iTempVar)

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

      !------------------------------------------------------------
      call HeatAD_reInitGridVars()
      !------------------------------------------------------------

      !------------------------------------------------------------
      call Simulation_adjustEvolution(dr_nstep, dr_dt, dr_simTime)
      !------------------------------------------------------------

      ! Fill GuardCells for temperature
      gcMask = .FALSE.
      gcMask(iTempVar) = .TRUE.
      call Grid_fillGuardCells(CENTER, ALLDIR, &
                               maskSize=NUNK_VARS + NDIM*NFACE_VARS, mask=gcMask)

      !------------------------------------------------------------
      call HeatAD_diffusion()
      call HeatAD_solve(dr_dt)
      !------------------------------------------------------------

      !------------------------------------------------------------
      call HeatAD_indicators()
      !------------------------------------------------------------

      !------------------------------------------------------------
      !- End Physics Sequence
      !------------------------------------------------------------
      call Grid_updateRefinement(dr_nstep, dr_simTime, gridChanged)

      !output a plotfile before the grid changes
      call Timers_start("IO_output")

      call IO_output(dr_simTime, &
                     dr_dt, dr_nstep + 1, dr_nbegin, endRunPl, PLOTFILE_AND_PARTICLEFILE)
      call Timers_stop("IO_output")

      if (gridChanged) dr_simGeneration = dr_simGeneration + 1

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

   return

end subroutine Driver_evolveAll
