!!****if* source/physics/Eos/EosMain/Gamma/eos_initGamma
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
!!  eos_initGamma
!!
!! SYNOPSIS
!!  
!!  call  eos_initGamma()
!!                 
!!
!! DESCRIPTION
!!
!!  This routine does ideal gamma law specific initialization
!!
!!  ARGUMENTS
!!
!!
!!***

#include "Eos.h"
subroutine eos_initGamma()

  use Eos_data, ONLY : eos_type
  use Eos_data, ONLY : eos_gamma
  use eos_idealGammaData, ONLY : eos_gammam1

  implicit none

  eos_type=EOS_GAM

  eos_gammam1 = 1.0/(eos_gamma-1.0)

  return
end subroutine eos_initGamma
