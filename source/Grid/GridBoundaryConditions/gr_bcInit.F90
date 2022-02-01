!!****if* source/Grid/GridBoundaryConditions/gr_bcInit
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
!!  gr_bcInit
!!
!! SYNOPSIS
!!
!!  gr_bcInit()
!!
!! DESCRIPTION
!!
!!  Initialize values for all data in the module gr_ptData,
!!  and allocate the scratch buffers
!!
!! ARGUMENTS
!!
!!
!!***

subroutine gr_bcInit()
  use gr_bcData, ONLY : gr_bcUseGravity, gr_bcEintSwitch

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

#include "Simulation.h"

  implicit none 

  gr_bcUseGravity = .false.
#ifdef GRAVITY
  call RuntimeParameters_get("useGravity",gr_bcUseGravity)
#endif

  gr_bcEintSwitch = 0.0
#ifdef FLASH_EOS
  call RuntimeParameters_get("eintSwitch",gr_bcEintSwitch)
#endif

  call gr_bcHseInit()

  return
end subroutine gr_bcInit
