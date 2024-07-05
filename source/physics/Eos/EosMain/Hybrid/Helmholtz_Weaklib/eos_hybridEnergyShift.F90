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
!! @brief eos_hybridEnergyShift implementation

!> @ingroup physics_sourceTerms_Burn
!!
!! @brief Compute the shift from internal to thermal energy
!!
!! @details
!!
!! This procedure will compute the total shift from internal
!! to thermal energy required for input into the Helmholtz
!! EOS.  This shift will account for the binding energy of each species
!! and the WeakLib offset relative to free neutrons
!!
!! @param ebin       Total binding energy in erg
!! @param massFrac   Mass fractions evolved by the simulation

subroutine eos_hybridEnergyShift(energyShift, vecLen, eosData, massFrac)
   use eos_hybridData, only: eos_hybDeltaE_WL, eos_hybBoverA, &
                             m_e, m_n, m_p, MeV2erg, N_A

#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

   implicit none

   integer, intent(in) :: vecLen
   real, dimension(vecLen), intent(out) :: energyShift
   real, dimension(vecLen,EOS_NUM),  intent(in) :: eosData
   real, dimension(NSPECIES,vecLen),  intent(in), optional :: massFrac

   integer :: i, k

   real, parameter :: conv = MeV2erg*N_A
   real, parameter :: delta = (m_n - m_p - m_e)*conv


   energyShift = eos_hybDeltaE_WL*conv

   ! Redundant yes, but I'd rather not have the conditional check
   ! in the loop over the vector
   if (present(massFrac) .and. (NSPECIES .gt. 0)) then
      ! Note: massFrac is indexed oppositely of eosData
      do k = 1, vecLen
         energyShift(k) = energyShift(k) - delta*eosData(k, EOS_YE) &
                          - sum(massFrac(:, k)*eos_hybBoverA)*conv
      end do ! k
   else
      do k = 1, vecLen
         energyShift(k) = energyShift(k) - delta*eosData(k, EOS_YE)
      end do ! k
   end if

end subroutine eos_hybridEnergyShift
