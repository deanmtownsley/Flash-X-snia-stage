!!****if* source/monitors/Profiler/ProfilerMain/mpihpm/Profiler_start
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
!!  Profiler_start
!!
!! SYNOPSIS
!!
!!  Profiler_start(character(IN) :: name, OR
!!                 integer(IN)   :: id)
!!                   
!!  
!! DESCRIPTION 
!!  
!!  Start profiling the section 'name' or 'id'
!!  
!! ARGUMENTS 
!!
!!  name - the name of the section to profile
!!  id - the integer id of the section to profile
!!
!!***

#include "Profiler.h"

subroutine Profiler_startName(name)
  use Profiler_data, ONLY : prf_evolutionName, prf_evolutionOnly
  use pr_interface, ONLY : pr_prof_control
  implicit none
  character (len=*), intent(in) :: name
  call hpm_start(name)
  if (trim(name) == prf_evolutionName .and. prf_evolutionOnly) then
     !summary_start zeroes out mpi_profile data collected up to this point.
     call summary_start()
     call pr_prof_control(PRF_START_PROFILER)
  end if
end subroutine Profiler_startName

subroutine Profiler_startId(id)
  use Driver_interface, ONLY : Driver_abort
  implicit none
  integer, intent(in) :: id
  call Driver_abort("Not yet implemented")
end subroutine Profiler_startId
