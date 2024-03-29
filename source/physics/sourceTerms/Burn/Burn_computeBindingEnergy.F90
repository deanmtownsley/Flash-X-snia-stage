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
!! @brief Burn_computeBindingEnergy stub

!> @ingroup physics_sourceTerms_Burn
!!
!! @brief Compute the total binding energy for species in the network
!!
!! @details
!! @anchor Burn_computeBindingEnergy_stub
!!
!! This function will compute the total binding energy offset to the fluid
!! internal energy.
!!
!! @param ebin       Total binding energy in erg
!! @param massFrac   Mass fractions evolved by the simulation
subroutine Burn_computeBindingEnergy(ebin, massFrac)
#include "Simulation.h"

   implicit none

   real, intent(out) :: ebin
   real, dimension(NSPECIES), intent(in) :: massFrac

   ebin = 0.0

   return
end subroutine Burn_computeBindingEnergy
