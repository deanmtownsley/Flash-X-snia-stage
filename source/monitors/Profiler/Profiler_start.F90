!!****f* source/monitors/Profiler/Profiler_start
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
implicit none
  character (len=*), intent(in) :: name

end subroutine Profiler_startName

subroutine Profiler_startId(id)
implicit none
  integer, intent(in) :: id

end subroutine Profiler_startId
