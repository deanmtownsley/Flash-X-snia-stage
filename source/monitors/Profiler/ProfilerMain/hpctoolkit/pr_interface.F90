!!****ih* source/monitors/Profiler/ProfilerMain/hpctoolkit/pr_interface
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
!!  pr_interface
!!
!! SYNOPSIS
!!  use pr_interface
!!
!! DESCRIPTION
!! This is the interface module for the profiler unit.
!!
!!***

module pr_interface
   interface
      subroutine pr_prof_control(mode) &
         bind(c, name='pr_prof_control')
         use iso_c_binding, ONLY: c_int
         integer(c_int), value, intent(in) :: mode
      end subroutine pr_prof_control
   end interface

   interface
      subroutine hpctoolkit_sampling_start() &
         bind(c, name='hpctoolkit_sampling_start')
      end subroutine hpctoolkit_sampling_start
   end interface

   interface
      subroutine hpctoolkit_sampling_stop() &
         bind(c, name='hpctoolkit_sampling_stop')
      end subroutine hpctoolkit_sampling_stop
   end interface
end module pr_interface
