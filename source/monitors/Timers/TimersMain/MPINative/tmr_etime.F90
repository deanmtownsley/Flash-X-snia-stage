!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_etime
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
!!   
!! tmr_etime
!!
!! SYNOPSIS
!!
!!  call tmr_etime(real(OUT) :: time)
!!
!! DESCRIPTION
!!
!!  Return the elapsed time.
!!
!! ARGUMENTS
!!
!!  time - elapsed time
!!
!!***
subroutine tmr_etime (time)
 
  implicit none 
  include 'mpif.h'
  real, intent(out) :: time

  time = MPI_WTime()
  return
end subroutine tmr_etime
