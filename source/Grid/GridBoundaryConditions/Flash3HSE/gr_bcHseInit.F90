!!****if* source/Grid/GridBoundaryConditions/Flash3HSE/gr_bcHseInit
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
!!  gr_bcHseInit
!!
!! SYNOPSIS
!!
!!  call gr_bcHseInit()
!!
!! DESCRIPTION
!!
!!  Initialization for a hydrostatic equilibrium (HSE) implementation
!!
!! ARGUMENTS
!!
!!
!!***

subroutine gr_bcHseInit()

  use gr_bcHseData, ONLY : gr_bcHseGravDirec, gr_bcHseDirection, gr_bcHseGravConst

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abort

#include "constants.h"
  implicit none 

  call RuntimeParameters_get('gconst', gr_bcHseGravConst)
  call RuntimeParameters_get("gdirec",gr_bcHseGravDirec)
  select case (gr_bcHseGravDirec)
  case('x','X')
     gr_bcHseDirection = IAXIS
  case('y','Y')
     gr_bcHseDirection = JAXIS
  case('z','Z')
     gr_bcHseDirection = KAXIS
  case default
     call Driver_abort('Runtime parameter "gdirec" only allows "x", "y", and "z".')
  end select

  return
end subroutine gr_bcHseInit
