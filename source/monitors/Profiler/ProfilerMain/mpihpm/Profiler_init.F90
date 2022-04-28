!!****if* source/monitors/Profiler/ProfilerMain/mpihpm/Profiler_init
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

#include "Profiler.h"

subroutine Profiler_init()
  use Profiler_data, ONLY : prf_evolutionOnly
  use pr_interface, ONLY : pr_prof_control
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  call RuntimeParameters_get("profileEvolutionOnly", prf_evolutionOnly)
  if (prf_evolutionOnly) then
     !If we wish to profile evolution and are using gprof then
     !we must disable gprof measurements first and then re-enable measurements
     !later when we are in the evolution section of FLASH.
     call pr_prof_control(PRF_DISABLE_PROFILER)
  end if
end subroutine Profiler_init
