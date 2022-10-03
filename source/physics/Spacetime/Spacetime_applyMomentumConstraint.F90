!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! @brief Spacetime_applyMomentumConstraint stub

!> @ingroup Spacetime
!!
!! @brief Apply the Hamiltonian constraint
!!
!! @details
!! @anchor Spacetime_applyMomentumConstraint_stub
!!
!! This procedure is responsible for applying the momentum constraint
!! if required by an implementation.
!!
!! @param  t  Time of the current solution of the evolved variables
subroutine Spacetime_applyMomentumConstraint(t)
   implicit none

   real, intent(in) :: t

   return
end subroutine Spacetime_applyMomentumConstraint
