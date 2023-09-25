!!****ih* source/physics/Eos/EosMain/hybrid/Helmholtz_WeakLib/eos_hybridInterface
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
!!     eos_hybridInterface
!!
!! SYNOPSIS
!!     use eos_hybridInterface
!!
!! DESCRIPTION
!!
!! This is an interface module for internal use of
!! the hybrid Eos implementation.
!!
!!***

module eos_hybridInterface

   implicit none

#include "Eos.h"
#include "Simulation.h"
#include "constants.h"

   interface
      subroutine eos_hybridSetFlag(mode, vecLen, eosData, massFrac, mask, vecBegin, vecEnd, eos_hybFlag)
         implicit none
         integer, INTENT(in) :: mode, vecLen
         real, INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData
         real, optional, INTENT(in), dimension(NSPECIES*vecLen) :: massFrac
         logical, optional, INTENT(in), target, dimension(EOS_VARS + 1:EOS_NUM) :: mask
         integer, optional, INTENT(in) :: vecBegin, vecEnd
         integer, optional, INTENT(out) :: eos_hybFlag(vecLen)
      end subroutine eos_hybridSetFlag
   end interface

end module eos_hybridInterface
