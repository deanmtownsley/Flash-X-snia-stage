!!****f* source/physics/ImBound/ImBound_init
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
!!  ImBound_init
!!
!!
!! SYNOPSIS
!!
!!  call ImBound_init(LOGICAL(IN) :: restart)
!!  
!! ARGUMENTS
!!
!!  restart - restart flag.
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which typically take values from runtime parameters.
!!  This must be called once by Driver_initAll.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine ImBound_init(restart)
  implicit none
  logical, INTENT(IN) :: restart
end subroutine ImBound_init
