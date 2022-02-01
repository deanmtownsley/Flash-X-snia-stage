!!****h* source/monitors/Profiler/Profiler_interface
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
!! This is the header file for the Profiler module that defines its
!! public interfaces.  Of course, there's no implementation of Profiler yet,
!! so calling these public interfaces doesn't do much good....
!!***

Module Profiler_interface

  interface Profiler_getSummary
     subroutine Profiler_getSummary(nIntervals)
       implicit none
       integer, intent(in) :: nIntervals
     end subroutine Profiler_getSummary
  end interface

  interface Profiler_init
     subroutine Profiler_init()
       implicit none
     end subroutine Profiler_init
  end interface

  interface Profiler_start
     subroutine Profiler_startName(name)
       implicit none
       character (len=*), intent(in) :: name
     end subroutine Profiler_startName

     subroutine Profiler_startId(id)
       implicit none
       integer, intent(in) :: id
     end subroutine Profiler_startId
  end interface

  interface Profiler_stop
     subroutine Profiler_stopName(name)
       implicit none
       character (len=*), intent(in) :: name
     end subroutine Profiler_stopName

     subroutine Profiler_stopId(id)
       implicit none
       integer, intent(in) :: id
     end subroutine Profiler_stopId
  end interface

end Module Profiler_interface
