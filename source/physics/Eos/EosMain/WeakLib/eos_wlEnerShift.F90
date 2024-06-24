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
!! @brief eos_wlEnerShift stub

!> @ingroup physics_eos_EosMain_WeakLib
!!
!! @brief Obtain the energy shift assumed by WeakLib
!!
!! @details
!! @anchor eos_wlEnerShift_stub
!!
!! This procedure can be used to obtain the energy shift that is assumed by WeakLib
!! It may be necessary to remove this shift when coupling to another Eos implementation,
!! e.g. when coupling to Helmholtz in the Hybrid Eos implementation, this energy shift
!! must be removed prior to passing the Eos data to Helmholtz
!!
!! @param energyShift   The energy shift to use
!!
subroutine eos_wlEnerShift(energyShift)
   implicit none

   real, intent(out) :: energyShift

   energyShift = 0.0

   return
end subroutine eos_wlEnerShift
