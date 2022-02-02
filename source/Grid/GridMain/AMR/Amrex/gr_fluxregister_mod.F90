!!****ih* source/Grid/GridMain/AMR/Amrex/gr_fluxregister_mod
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
!!  gr_fluxregister_mod
!!
!! SYNOPSIS
!!
!!  use gr_fluxregister_mod,ONLY: gr_fluxregister_t, gr_fluxregisterBuild, gr_fluxregisterDestroy
!!
!! DESCRIPTION
!!
!!  This module exists to make available the public type dieclarations
!!  and subroutines from either amrex_flash_fluxregister_module or
!!  amrex_fluxregister_module under unified names.
!!
!!  Which of the two is used defends on the preprocessor symbol
!!  USE_AMREX_FLASHFLUXREGISTER .
!!
!! NOTES
!!
!!  This module is only useful for Grid implementations that
!!  implement SPFS (semi-permanent flux storage) in a form of "flux
!!  registers", i.e., Amrex.
!!***

#include "Simulation.h"

module gr_fluxregister_mod
#ifdef USE_AMREX_FLASHFLUXREGISTER
  use amrex_flash_fluxregister_module, ONLY: gr_fluxregister_t      => amrex_flash_fluxregister
  use amrex_flash_fluxregister_module, ONLY: gr_fluxregisterBuild   => amrex_flash_fluxregister_build
  use amrex_flash_fluxregister_module, ONLY: gr_fluxregisterDestroy => amrex_flash_fluxregister_destroy
#else
  use amrex_fluxregister_module, ONLY: gr_fluxregister_t      => amrex_fluxregister
  use amrex_fluxregister_module, ONLY: gr_fluxregisterBuild   => amrex_fluxregister_build
  use amrex_fluxregister_module, ONLY: gr_fluxregisterDestroy => amrex_fluxregister_destroy

#endif

  implicit none

end module gr_fluxregister_mod
