!!****if* source/monitors/Profiler/ProfilerMain/hpctoolkit/Profiler_stop
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
!!  Profiler_stop
!!
!! SYNOPSIS
!!
!!  Profiler_stop( character(IN) :: name, OR
!!                 integer(IN)   :: id)
!!
!!
!! DESCRIPTION
!!
!!  Stop profiling the section 'name' or 'id'
!!
!! ARGUMENTS
!!
!!  name - the name of the section to profile
!!  id - the integer id of the section to profile
!!
!!***

subroutine Profiler_stopName(name)
   use Profiler_data, ONLY: prf_profilerIsOn, prf_groupName
   use pr_interface, ONLY: hpctoolkit_sampling_stop
   use Driver_interface, ONLY: Driver_abort
   implicit none
   character(len=*), intent(in) :: name
   character(len=200) :: errorMessage
   if (trim(prf_groupName) == name) then
      if (prf_profilerIsOn) then
         call hpctoolkit_sampling_stop()
         prf_profilerIsOn = .FALSE.
      else
         write (errorMessage, *) "[Profiler_stop] Cannot match ", trim(name), " with a previously started Profiler"
         call Driver_abort(trim(errorMessage))
      end if
   end if
end subroutine Profiler_stopName

subroutine Profiler_stopId(id)
   use Driver_interface, ONLY: Driver_abort
   implicit none
   integer, intent(in) :: id
   call Driver_abort("[Profiler_stop] Not yet implemented using integer ID as argument")
end subroutine Profiler_stopId
