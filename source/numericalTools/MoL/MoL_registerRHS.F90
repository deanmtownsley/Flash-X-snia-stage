!> @file source/numericalTools/MoL/MoL_registerRHS.F90
!!
!! @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! @brief MoL_registerRHS stub
!! @ingroup MoL

!> @brief Register a procedure responsible for calculating RHS terms
!!
!! @param rhsType  RHS-type identifier
!! @param rhsFunc  Procedure that will calculate the RHS terms
!! @anchor MoL_registerRHS_stub
!!
!! @pre `rhsType` is a valid MoL RHS identifier as defined in @ref MoL.h
!! @pre `rhsFunc` is a callable procedure of the current time: `call rhsFunc(t)`
!!
!! @sideffect MoL's procedure pointer for the specified RHS-type will be set to the
!!            provided procedure
!!
!! @returns None
!!
!! @details
!!    Valid RHS types include (defined in MoL.h):
!!       - `MOL_RHS_EXPLICIT`  : RHS for (slow) explicit terms
!!       - `MOL_RHS_IMPLICIT`  : RHS for (slow) implicit terms
!!       - `MOL_RHS_FAST`      : RHS for (fast) explicit terms
!!
!! @warning Will trigger Flash-X to abort if an invalid RHS-type is specified
!!
!! @ingroup MoL
subroutine MoL_registerRHS(rhsType, rhsFunc)
   use MoL_functionTypes, only: MoL_rhs_t

   implicit none

   integer, intent(in) :: rhsType
   procedure(MoL_rhs_t) :: rhsFunc

   return
end subroutine MoL_registerRHS
