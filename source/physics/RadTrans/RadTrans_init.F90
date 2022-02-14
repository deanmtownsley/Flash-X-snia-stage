!!****f* source/physics/RadTrans/RadTrans_init
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
!!  NAME 
!!
!!  RadTrans_init
!!
!!  SYNOPSIS
!!
!!  call RadTrans_init()
!!
!!  DESCRIPTION 
!!    Initialize radiative transfer unit
!!
!! ARGUMENTS
!!
!!
!!***
subroutine RadTrans_init()

  ! Stub implementation

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abort

  implicit none

  logical, save :: testUseRadTrans

  !! It is a failure to invoke the stub when useRadTrans is set TRUE.

  call RuntimeParameters_get ("useRadTrans", testUseRadTrans)
  if (testUseRadTrans) then
     call Driver_abort("RadTrans unit seems not to be compiled in, and the RadTrans_init stub does not &
          &allow the value of useRadTrans to be TRUE.")
  end if

  return
end subroutine RadTrans_init
