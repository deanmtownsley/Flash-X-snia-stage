!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!> @ingroup physics_Eos
!!
!! @brief This is an interface module for internal use of
!!        the hybrid Eos implementation.
!!

module eos_hybridInterface

   implicit none

#include "Eos.h"
#include "Simulation.h"
#include "constants.h"

   interface
      subroutine eos_hybridSetFlag(mode, vecLen, eosData, vecBegin, vecEnd, eos_hybFlag)
         implicit none
         integer, intent(in) :: mode, vecLen
         real, intent(in), dimension(vecLen,EOS_VARS) :: eosData
         integer, intent(in) :: vecBegin, vecEnd
         integer, intent(out) :: eos_hybFlag(vecLen)
      end subroutine eos_hybridSetFlag
   end interface

   interface
      subroutine eos_hybridEnergyShift(energyShift, vecLen, eosData, massFrac)
         implicit none

         integer, intent(in) :: vecLen
         real, dimension(vecLen), intent(out) :: energyShift
         real, dimension(vecLen,EOS_VARS), target, intent(in) :: eosData
         real, dimension(NSPECIES,vecLen), target, intent(in), optional :: massFrac
      end subroutine eos_hybridEnergyShift
   end interface

end module eos_hybridInterface
