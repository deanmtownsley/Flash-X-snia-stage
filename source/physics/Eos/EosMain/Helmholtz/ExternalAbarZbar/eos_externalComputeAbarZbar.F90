!!****if* source/physics/Eos/EosMain/Helmholtz/ExternalAbarZbar/eos_externalComputeAbarZbar
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
!!  eos_externalComputeAbarZbar
!!
!! SYNOPSIS
!!
!!  call eos_externalComputeAbarZbar(real(in), dimension(:,:)  :: solnscalars,
!!                                   real(out), dimension(:)  :: abardata,
!!                                   real(out), dimension(:)  :: zbardata)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!! Klaus Weide   2013
!!
!!  This routine private to the Eos unit serves to implement callbacks
!!  for the Eos/Helmholtz/ExternelAbarZbar EOS implementation.
!!  Code units that implement ways for computing Abar and Zbar from
!!  the solutions state (or otherwise) are polled here by calling their
!!  public UNITNAME_computeAbarZbar interfaces.
!!
!!  It is assumed that not more than one polled units have implementations
!!  included in the simulation that actually provide Abar and Zbar.
!!  If several units do provide the data, the last one polled here will win.
!!
!! ARGUMENTS
!!
!!   solnscalars : scalars of the solution 
!!
!!   abardata : abar info
!!
!!   zbardata : zbar info
!!
!!
!!
!!***

#include "Simulation.h"
subroutine eos_externalComputeAbarZbar(solnScalars, abarData, zbarData)

  use Flame_interface, ONLY: Flame_computeAbarZbar
  use Burn_interface,  ONLY: Burn_computeAbarZbar

  implicit none

  real, intent(in),  dimension(:,:)  :: solnScalars
  real, intent(out), dimension(:)    :: abarData, zbarData

  call Flame_computeAbarZbar(solnScalars, abarData, zbarData)

#ifdef FLASH_SOURCEBURN
  call Burn_computeAbarZbar(solnScalars, abarData, zbarData)
#endif

end subroutine eos_externalComputeAbarZbar
