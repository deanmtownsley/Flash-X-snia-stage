!!****if* source/monitors/Profiler/ProfilerMain/Profiler_init
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
!!  Profiler_init
!!
!! SYNOPSIS
!!
!!  Profiler_init()
!!
!!
!! DESCRIPTION
!!
!!  Initialize the profiler unit
!!
!! ARGUMENTS
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine Profiler_init()
   use Profiler_data, ONLY: prf_profilerIsOn, prf_groupName, prf_evolutionOnly, prf_meshMe
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use Driver_interface, ONLY: Driver_abort, Driver_getMype
   use Logfile_interface, ONLY: Logfile_stamp

   implicit none

   call Driver_getMype(MESH_COMM, prf_meshMe)

   ! Initialize profiler specific data and read runtime parameters.
   prf_profilerIsOn = .FALSE.
   call RuntimeParameters_get("profileEvolutionOnly", prf_evolutionOnly)
   call RuntimeParameters_get("profileGroupName", prf_groupName)

   ! Perform a sanity check on runtime parameters and abort if
   ! inconsistencies are found in their definitions
   if (.not. prf_evolutionOnly .and. trim(prf_groupName) == "FLASHX_EVOLUTION") then
      call Driver_abort("[Profiler_init] profileGroupName is FLASHX_EVOLUTION but profileEvolutionOnly is .FLASE.")
   else if (prf_evolutionOnly .and. trim(prf_groupName) /= "FLASHX_EVOLUTION") then
      call Driver_abort("[Profiler_init] profileEvolutionOnly is .TRUE. but profileGroupName is not FLASHX_EVOLUTION")
   end if

   call Logfile_stamp(trim(prf_groupName), "[Profiler_init]")

   if (prf_meshMe .eq. MASTER_PE) then
      print *, "Profiling group name: ", trim(prf_groupName)
   end if

end subroutine Profiler_init
