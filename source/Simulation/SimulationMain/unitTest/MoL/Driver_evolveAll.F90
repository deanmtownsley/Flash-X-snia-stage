!!****if* source/Simulation/SimulationMain/unitTest/MoL/Driver_evolveAll
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
!! MoL compatible version of Driver_evolveAll
!!
!! NOTES
!!
!!***
subroutine Driver_evolveAll()

   use Driver_data, ONLY: dr_globalMe, dr_nbegin, &
                          dr_nend, dr_dt, &
                          dr_tmax, dr_simTime, &
                          dr_nstep, dr_dtOld, dr_dtNew, &
                          dr_simGeneration
   use Driver_interface, ONLY: Driver_computeDt, Driver_driftUnk
   use Logfile_interface, ONLY: Logfile_stamp, Logfile_close
   use Timers_interface, ONLY: Timers_start, Timers_stop, &
                               Timers_getSummary
   use Profiler_interface, ONLY: Profiler_start, Profiler_stop
   use Grid_interface, ONLY: Grid_updateRefinement

   use MoL_interface, only: MoL_regrid, &
                            MoL_advance, &
                            MoL_releaseFunctions

   use dr_molInterface, only: dr_molRegisterFunctions, &
                              dr_molPreEvolve, &
                              dr_molPostTimeStep, &
                              dr_molPostRegrid

   use sim_molInterface, only: sim_molVerifySolution

#include "Simulation.h"
#include "constants.h"

   implicit none

   ! for logfile output
   character(len=MAX_STRING_LENGTH), dimension(4, 2) :: strBuff
   character(len=15) :: numToStr

   logical :: gridChanged
   logical :: endRunWallClock !Should we end our run on this iteration, based on wall clock time?
   logical :: shortenedDt !Is the last timestep being shortened to reach dr_tmax?

   integer, parameter :: driftUnk_flags = DRIFT_NO_PARENTS

   logical :: valid
   real :: maxError

   character(len=20) :: fileName
   integer, parameter        :: fileUnit = 2
   integer, dimension(4) :: prNum
   integer :: temp, i

   temp = dr_globalMe

   do i = 1, 4
      prNum(i) = mod(temp, 10)
      temp = temp/10
   end do
   filename = "unitTest_"//char(48 + prNum(4))//char(48 + prNum(3))// &
              char(48 + prNum(2))//char(48 + prNum(1))

   open (fileUnit, file=fileName)
   write (fileUnit, '("P",I0)') dr_globalMe

  !! MoL needs a regrid call to setup initial storage for intermediate states
  !! This must occur after ALL *_init calls that have MoL_registerVariable calls
   call MoL_regrid

   call Logfile_stamp('Entering evolution loop', '[Driver_evolveAll]')
   call Profiler_start("FLASH_evolution")
   call Timers_start("evolution")

   ! Not necessarily the ideal location for this, but it needs to be done before the
   ! the evolution loop begins, but after all initialization is complete
   call dr_molPreEvolve(dr_simTime)

   do dr_nstep = dr_nBegin, dr_nend

      call dr_shortenLastDt(dr_dt, dr_simTime, dr_tmax, shortenedDt, 1)
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

         call Logfile_stamp(strBuff(1:3, :), 3, 2, "step")

      end if

      call Driver_driftUnk(__FILE__, __LINE__, driftUnk_flags)

     !! ================ !!
     !!   MoL evolution  !!
     !! ================ !!

      ! Each MoL-based simulation will register RHS, postUpdate, etc. calls w/ MoL
      call dr_molRegisterFunctions

      ! Take a single time step from t to t+dt
      call MoL_advance(dr_simTime, dr_dt)

      ! Reset registered functions so other physics units can utilitize MoL below
      call MoL_releaseFunctions

      ! This has to occur after MoL advances - MoL requires the time at the start
      ! of a timestep, not the end...
      dr_simTime = dr_simTime + dr_dt

     !! This will notify MoL-evolved physics units that timestep is complete
      call dr_molPostTimeStep(dr_simTime)

      dr_dtOld = dr_dt

      dr_simGeneration = 0

      !----
      !- End Physics Sequence
      !--------------------------------------------------------------------

      call Timers_start("Grid_updateRefinement")
      call Grid_updateRefinement(dr_nstep, dr_simTime, gridChanged)
      call Timers_stop("Grid_updateRefinement")

      ! Update MoL intermediate state storage if the grid changed
      ! Also notify physics units - if they track both conserved
      ! and primitive variables in UNK, the refined primitives will
      ! not be consistent with a con2prim call on the refined
      ! conserved variables.  A con2prim call may be necessary
      ! here, as the assumed strictly-Newtonian expressions for
      ! con2prim and prim2con in Grid are incorrect for any
      ! relativistic formulation
      if (gridChanged) then
         call MoL_regrid
         call dr_molPostRegrid(dr_simTime)
      end if

      if (gridChanged) dr_simGeneration = dr_simGeneration + 1

      ! calculate new
      call Timers_start("Driver_computeDt")
      call Driver_computeDt(dr_nbegin, dr_nstep, &
                            dr_simTime, dr_dtOld, dr_dtNew)
      call Timers_stop("Driver_computeDt")

      ! store new
      dr_dt = dr_dtNew

     !!*****************************************************************************
     !!  Evolution Loop -- check termination conditions
     !!*****************************************************************************

      !Exit if this step was handled specially as the last step
      if (shortenedDt) exit

     !! the simulation ends before nend iterations if
     !!  (i)   the simulation time is greater than the maximum time (tmax)
     !!  (ii)  the redshift falls below the minimum redshift
     !!        (also called redshiftFinal)
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

   call sim_molVerifySolution(dr_simTime, dr_dt, valid, maxError)
   if (dr_globalMe == MASTER_PE) then
      print *, "MoL unit test passed?", valid
      print *, "Max error: ", maxError
   end if

   if (valid) then
      write (fileUnit, '(A)') 'SUCCESS all results conformed with expected values.'
   else
      write (fileUnit, '(A)') 'FAILURE'
   end if

   close (fileUnit)

   !The value of dr_nstep after the loop is (dr_nend + 1) if the loop iterated for
   !the maximum number of times.  However, we need to retain the value that
   !dr_nstep had during the last loop iteration, otherwise the number for nstep
   !that will be stored in a final checkpoint file will be wrong.
   dr_nstep = min(dr_nstep, dr_nend)

  !!******************************************************************************
  !! End of Evolution Loop
  !!******************************************************************************

   call Timers_stop("evolution")
   call Profiler_stop("FLASH_evolution")
   call Logfile_stamp('Exiting evolution loop', '[Driver_evolveAll]')
   call Timers_getSummary(max(0, dr_nstep - dr_nbegin + 1))
   call Logfile_stamp("FLASH run complete.", "LOGFILE_END")
   call Logfile_close()

   return

end subroutine Driver_evolveAll
