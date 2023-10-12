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
      subroutine eos_hybridSetFlag(mode, vecLen, eosData, vecBegin, vecEnd, eos_hybFlag)
         implicit none
         integer, intent(in) :: mode, vecLen
         real, intent(in), dimension(EOS_NUM*vecLen) :: eosData
         integer, intent(in) :: vecBegin, vecEnd
         integer, intent(out) :: eos_hybFlag(vecLen)
      end subroutine eos_hybridSetFlag
   end interface

   interface
      subroutine eos_hybridEnergyShift(energyShift, vecLen, eosData, massFrac)
         implicit none

         integer, intent(in) :: vecLen
         real, dimension(vecLen), intent(out) :: energyShift
         real, dimension(vecLen*EOS_NUM), target, intent(in) :: eosData
         real, dimension(vecLen*EOS_NUM), target, intent(in), optional :: massFrac
      end subroutine eos_hybridEnergyShift
   end interface

   interface
      subroutine eos_hybridHelmholtz(mode, eosData, mask, massFrac)
         implicit none
         integer, intent(in) :: mode
         real, dimension(EOS_NUM), target, intent(inout) :: eosData
         logical, dimension(EOS_VARS + 1:EOS_NUM), target, intent(in) :: mask
         real, dimension(NSPECIES), target, intent(in), optional :: massFrac
      end subroutine eos_hybridHelmholtz
   end interface

   interface
      subroutine eos_hybridWeakLib(mode, eosData, mask, massFrac)
         implicit none
         integer, intent(in) :: mode
         real, dimension(EOS_NUM), target, intent(inout) :: eosData
         logical, dimension(EOS_VARS + 1:EOS_NUM), target, intent(in) :: mask
         real, dimension(NSPECIES), target, intent(in), optional :: massFrac
      end subroutine eos_hybridWeakLib
   end interface

end module eos_hybridInterface
