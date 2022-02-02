!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_create
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
!!  tmr_create
!!
!! SYNOPSIS
!!
!!  tmr_create(character(IN) :: name, 
!!                 integer(OUT), optional :: i)
!!
!! DESCRIPTION
!!  Create a new timer, and optionally return the
!!  integer key to that timer
!!
!! ARGUMENTS
!!  name --    a string containing the name of the section to time
!!  i --       an optional integer parameter that will return the
!!              integer key to the timer.
!!
!! PARAMETERS
!!
!!***

subroutine tmr_create (name,i)

  use Timers_data, ONLY: tmr_timerInvalid
  implicit none

! Arguments
  character(len=*), intent(IN)  :: name
  integer,optional, intent(OUT) :: i
! Local variables
  integer          :: id
  
#ifdef NOOP
  return
#endif
  
  call tmr_findTimerIndex (name,.FALSE., id)
  
  if (id == tmr_timerInvalid) then
     call tmr_findTimerIndex (name,.TRUE., id)
  endif
  
  if (present(i)) i = id
  
  return
  
end subroutine tmr_create
