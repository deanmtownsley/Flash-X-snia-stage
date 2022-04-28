!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_lookupIndex
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
!!   tmr_lookupIndex 
!!
!! SYNOPSIS
!!
!!   call tmr_lookupIndex( character(len=*), intent(IN)  :: name, 
!!                              int, intent(OUT)         :: index)
!!
!! DESCRIPTION
!!   find the integer key (for faster use) for a given name
!!
!! ARGUMENTS
!!
!!   name --  a string containing the name of a timer
!!   index -- the returned timer index
!!
!! PARAMETERS
!!
!!***

subroutine tmr_lookupIndex(name, index)

  implicit none  

  character(len=*), intent(IN) :: name
  integer, intent(out) :: index
  
  call tmr_findTimerIndex(name,.FALSE., index)
  
  return
  
end subroutine tmr_lookupIndex
