!!****f* source/monitors/Timers/Timers_stop
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
!!   Timers_stopIndex - stop a timer given an integer key
!!
!! SYNOPSIS
!!
!!   Timers_stopIndex(integer(IN) :: i)
!!
!! DESCRIPTION
!!   Stop timing a timer specified by a supplied integer key
!!
!! ARGUMENTS
!!   i --   an integer key specifiying the timer to stop
!!
!! PARAMETERS
!!
!!***

recursive subroutine Timers_stopIndex (i)

  implicit none  
  integer, intent(in) :: i
  return
end subroutine Timers_stopIndex

!!****if* source/monitors/Timers/Timers_stopString
!!
!! NAME
!!
!!   Timers_stopString
!!
!! SYNOPSIS
!!   Timers_stopString(character(IN) :: name(len=*))
!!
!! DESCRIPTION
!!   Stop timing a timer specified by a supplied name
!!
!! ARGUMENTS
!!   name --   a name specifiying the timer to stop
!!
!! PARAMETERS
!!
!!***
subroutine Timers_stopString(name)

  implicit none  

  character(len=*), intent(IN)   :: name
  return
end subroutine Timers_stopString
