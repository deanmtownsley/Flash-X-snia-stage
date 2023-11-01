!!****if* source/monitors/Profiler/ProfilerMain/hpctoolkit/Profiler_start
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

subroutine Profiler_startName(name)
   use Profiler_data, ONLY: prf_profilerIsOn, prf_groupName
   use pr_interface, ONLY: hpctoolkit_sampling_start
   use Driver_interface, ONLY: Driver_abort
   implicit none
   character(len=*), intent(in) :: name
   character(len=200) :: errorMessage
   if (trim(prf_groupName) == name) then
      if (.not. prf_profilerIsOn) then
         prf_profilerIsOn = .TRUE.
         call hpctoolkit_sampling_start()
      else
         write (errorMessage, *) "[Profiler_start] Cannot start Profiler for ", trim(name), &
            " because a previous Profiler is already ON"
         call Driver_abort(trim(errorMessage))
      end if
   end if
end subroutine Profiler_startName

subroutine Profiler_startId(id)
   use Driver_interface, ONLY: Driver_abort
   implicit none
   integer, intent(in) :: id
   call Driver_abort("[Profiler_start] Not yet implemented using integer ID as argument")
end subroutine Profiler_startId
