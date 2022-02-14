!!****f* source/monitors/Timers/Timers_init
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
!!  Timers_init
!!
!! SYNOPSIS
!!
!!  Timers_init(real(OUT) :: initialWCTime)
!!  
!! DESCRIPTION 
!!  
!!  Initialize the timer data structures.  This will
!!  essentially delete all information previously gathered by all timers
!!  and make it safe to start timers from scratch.  In the middle of
!!  a run, for instance, this could be called once per timestep along with
!!  Timers_getSummary to get timer summary information for each timestep.
!!  
!! ARGUMENTS 
!!
!!  initialWCTime -- the initial wall clock time when this was called. 
!! 
!!***

subroutine Timers_init( initialWCTime)

  implicit none
  include "mpif.h"
  real, intent(out) :: initialWCTime

  initialWCTime = MPI_WTime()

end subroutine Timers_init
