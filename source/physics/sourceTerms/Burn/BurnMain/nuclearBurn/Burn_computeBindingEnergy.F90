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
!! @brief Burn_computeBindingEnergy implementation

!> @ingroup physics_sourceTerms_Burn
!!
!! @brief Implements Burn_computeBindingEnergy
!!
!! @details
!! @stubref{Burn_computeBindingEnergy}
!!
!! Based on BANG version in networkPlusHybrid branch and misc. course notes
!! from MSU's PHY 981/2/3
subroutine Burn_computeBindingEnergy(ebin, massFrac)
   use Burn_data, only: xmass, ymass, bion, zion, aion
   use Burn_dataEOS, only: bye, abar, ytot1
   use bn_interface, only: bn_azbar

#include "Simulation.h"

   implicit none

   real, intent(out) :: ebin
   real, dimension(NSPECIES), intent(in) :: massFrac

   ! From CODATA 2018
   real, parameter :: m_u = 931.49410242 ! MeV
   real, parameter :: m_n = 939.56542052 ! MeV
   real, parameter :: m_p = 938.27208816 ! MeV
   real, parameter :: m_e = 0.51099895   ! MeV

   real, parameter :: delta_n = m_n - m_u ! MeV
   real, parameter :: delta_p = m_p - m_u ! MeV

   real, parameter :: N_A = 6.02214076e23 ! mol^-1

   real, parameter :: MeV2erg = 1.602176634e-6 ! erg / MeV

   real, parameter :: conv = MeV2erg*N_A ! erg / mol

   integer :: s

   do s = 1, NSPECIES
      xmass(s) = massFrac(s)
   end do ! s

   call bn_azbar()

   ebin = 0.0

   do s = 1, NSPECIES
      ebin = ebin + ymass(s)*(bion(s) - zion(s)*delta_p - (aion(s) - zion(s))*delta_n)*conv
   end do ! s

   ! TODO: conditional check on whether the network tracks neutrons as a species (e.g. XNet)
   ebin = ebin + (bye*m_e + ytot1*delta_n)*conv
end subroutine Burn_computeBindingEnergy
