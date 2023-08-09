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
   use Driver_interface, only: Driver_getMype
   use Timers_interface, only: Timers_getSummary, Timers_start, Timers_stop
   use Logfile_interface, only: Logfile_stamp

   use Spacetime_interface, only: Spacetime_unitTest

   use iso_fortran_env, only: stdout => output_unit

#include "constants.h"

   implicit none

   logical ::  perfect = .true.

   integer :: procID

   call Driver_getMype(GLOBAL_COMM, procID)

   write (stdout, "(a,i0.4,a)") "Proc #", procID, " - Spacetime unit test started"

   perfect = .true.

   call Timers_start("testing")
   call Spacetime_unitTest(stdout, perfect)
   call Timers_stop("testing")

   if (perfect) then
      write (stdout, "(a,i0.4,a)") "Proc #", procID, " - Spacetime unit test PASSED"
   else
      write (stdout, "(a,i0.4,a)") "Proc #", procID, " - Spacetime unit test FAILED"
   end if

   call Timers_getSummary(0)

   call Logfile_stamp("Flash-X Spacetime unit test run complete.", "LOGFILE_END")
end subroutine Driver_evolveAll
