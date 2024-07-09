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
!! @brief  This is the data module for the hybrid (Helmholtz+WeakLib) Eos implementation.
!!         It stores all the runtime parameters, and all the unit scope
!!         data. Some of the unit scope data is fecthed by the wrapper layer
!!         from elsewhere in the code and some is local unit data common
!!         multiple functions in the unit
!!

module eos_hybridData
#include "Simulation.h"

   implicit none

   real, save :: eos_hybTransitionDensHi
   real, save :: eos_hybTransitionDensLo

   ! These were in old hybridData
   real, parameter :: amu_cgs = 1.66053873e-24
   real, parameter :: kb_erg = 1.380658e-16

   real, parameter :: ergPerKg_to_kbPerBaryon = amu_cgs/kb_erg
   real, parameter :: kbPerBaryon_to_ergPerKg = 1./ergPerKg_to_kbPerBaryon

   ! Need these onces from various Units
   real, save :: eos_hybDeltaE_WL
   real, dimension(NSPECIES), save :: eos_hybBoverA

   ! Ideally, these should be set from the values in the PhysicalConstants unit,
   ! but not all of these are currently available.  For now, these will be set
   ! to the CODATA 2018 values, but varying CODATA sets are present throughout
   ! other portions of the code, so these will need to be reconciled with those
   ! used in XNet, WeakLib, etc.

   ! Masses
   real, parameter :: m_u = 931.49410242 ! MeV /  c^2
   real, parameter :: m_n = 939.56542052 ! MeV /  c^2
   real, parameter :: m_p = 938.27208816 ! MeV /  c^2
   real, parameter :: m_e = 0.51099895   ! MeV /  c^2

   real, parameter :: delta_n = m_n - m_u       ! MeV / c^2
   real, parameter :: delta_p = m_p + m_e - m_u ! MeV / c^2

   real, parameter :: N_A = 6.02214076e23 ! mol^-1

   real, parameter :: MeV2erg = 1.602176634e-6  ! erg / MeV
   real, parameter :: erg2MeV = 1.0/MeV2erg     ! MeV / erg
end module eos_hybridData
