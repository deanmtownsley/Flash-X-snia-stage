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
!! @brief MoL_registerRHS stub

!> @ingroup MoL
!!
!! @brief Register a procedure responsible for calculating RHS terms
!!
!! @details
!! @anchor MoL_registerRHS_stub
!!
!! Valid RHS types include (defined in MoL.h):
!!    - `MOL_RHS_EXPLICIT`  : RHS for (slow) explicit terms
!!    - `MOL_RHS_IMPLICIT`  : RHS for (slow) implicit terms
!!    - `MOL_RHS_FAST`      : RHS for (fast) explicit terms
!!
!! @param rhsType  RHS-type identifier
!! @param rhsFunc  Procedure that will calculate the RHS terms
subroutine MoL_registerRHS(rhsType, rhsFunc)
   use MoL_functionTypes, only: MoL_rhs_t

   implicit none

   integer, intent(in) :: rhsType
   procedure(MoL_rhs_t) :: rhsFunc

   return
end subroutine MoL_registerRHS
