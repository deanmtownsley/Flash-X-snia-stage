!!****if* source/physics/Eos/EosMain/hybrid/Helmholtz_WeakLib/eos_hybridEnergyShift
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  eos_hybridEnergyShift
!!
!! SYNOPSIS
!!
!!  call eos_hybridEnergyShift(integer(IN)  :: mode,
!!                         integer(IN)  :: vecLen,
!!                         real(INOUT)  :: eosData(vecLen*EOS_NUM),
!!               optional, real(IN)     :: massFrac(vecLen*NSPECIES),
!!               optional, logical(IN)  :: mask(EOS_VARS+1:EOS_NUM))
!!
!! DESCRIPTION
!!
!!  This procedure prepares and calls the Helmholtz Eos
!!
!! ARGUMENTS
!!
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points for each input variable
!!
!!  eosData  : This array is the data structure through which variable values are
!!             passed in and out of the Eos routine. The arrays is sized as
!!             EOS_NUM*vecLen. EOS_NUM, and individual input and output
!!             Eos variables are defined in Eos.h. The array is organizes such that
!!             the first 1:vecLen entries represent the first Eos variable, vecLen+1:
!!             2*vecLen represent the second Eos variable and so on.
!!
!!  massFrac : If present, contains the mass fractions of the species included in
!!             the simulation. The array is sized as NSPECIES*vecLen.
!!             If you are working with electron abundance mass scalars, then you
!!             do not necessarily have to have mass fractions.
!!             Thus, for the Ye formulation this value is never used!
!!
!!  mask     : Mask is a logical array the size of EOS_DERIVS (number
!!              of partial derivatives that can be computed, defined in
!!              Eos.h), where each index represents a specific partial derivative
!!              that can be calculated by the Eos unit. A .true. value in mask
!!              results in the corresponding derivative being calculated and
!!              returned. It should preferably be dimensioned as
!!              mask(EOS_VARS+1:EOS_NUM) in the calling routine
!!              to exactly match the arguments declaration in Eos Unit.
!!             Note that the indexing of mask does not begin at 1, but rather at one past
!!             the number of variables.
!!
!!***
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
   real, dimension(vecLen*EOS_NUM), target, intent(in) :: eosData
   real, dimension(NSPECIES*vecLen), target, intent(in), optional :: massFrac

   real, dimension(:, :), pointer :: eosData_ptr, massFrac_ptr
   integer :: i, k

   real, parameter :: conv = MeV2erg*N_A
   real, parameter :: delta = (m_n - m_p - m_e)*conv

   eosData_ptr(1:vecLen, 1:EOS_NUM) => eosData

   energyShift = eos_hybDeltaE_WL*conv

   ! Redudant yes, but I'd rather not have the conditional check
   ! in the loop over the vector
   if (present(massFrac) .and. (NSPECIES .gt. 0)) then
      ! Note: massFrac is indexed oppositely of eosData
      massFrac_ptr(1:NSPECIES, 1:vecLen) => massFrac

      do k = 1, vecLen
         energyShift(k) = energyShift(k) + delta*eosData_ptr(k, EOS_YE) &
                          + sum(massFrac_ptr(:, k)*eos_hybBoverA)*conv
      end do ! k
   else
      do k = 1, vecLen
         energyShift(k) = energyShift(k) + delta*eosData_ptr(k, EOS_YE)
      end do ! k
   end if

end subroutine eos_hybridEnergyShift
