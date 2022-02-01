!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_getMaxCallStackDepth
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
!!   tmr_getMaxCallStackDepth
!!
!! SYNOPSIS
!!   tmr_getMaxCallStackDepth(integer(OUT) :: value)
!!
!! DESCRIPTION
!!  Accessor function for how deep the timers can nest at runtime, which is 
!!  in general, a fixed depth.
!!
!! ARGUMENTS
!!   value -- the max depth for nesting timers.
!!
!! PARAMETERS
!!
!!***

subroutine tmr_getMaxCallStackDepth(value)

  use Timers_data, ONLY: tmr_maxCallStackDepth

  implicit none

  integer, intent(OUT) :: value

  value = tmr_maxCallStackDepth
  return

end subroutine tmr_getMaxCallStackDepth
