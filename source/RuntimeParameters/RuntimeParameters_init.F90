!!****f* source/RuntimeParameters/RuntimeParameters_init
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
!!  RuntimeParameters_init
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_init(logical(out) :: restart)
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the data need in the Runtime Parameters
!!  Unit.  Reads the parameter file, usually flash.par,
!!  and broadcasts parameters to the other processors.
!!    
!!
!! ARGUMENTS
!!
!! restart - true if run is restarted from checkpoint, false if starting
!!           from scratch
!!           
!!
!!
!!***



subroutine RuntimeParameters_init( restart)

implicit none

  logical, intent(out) :: restart

  restart = .false.

end subroutine RuntimeParameters_init
