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
!! @brief ml_calcRHS stub

!> @ingroup MoLPrivate
!!
!! @brief Calculate specified RHS type and store in requested struct
!!
!! @details
!! @anchor ml_calcRHS_stub
!!
!! This procedure will properly account for fallthrough cases
!! for each integrator, e.g. an explicit-only integrator will
!! call all available registered RHS procedures to add into
!! the RHS terms used during each integration stage.
!!
!! Valid MoL data structures include (defined in @ref MoL.h):
!! - `MOL_RHS_EXPLICIT`
!! - `MOL_RHS_IMPLICIT`
!! - `MOL_RHS_FAST`
!!
!! @param rhsType    The type of RHS
!! @param rhsStruct  MoL memory data-struct to store RHS in
!! @param t          The time of the RHS is to be evaluated at
!! @param dtWeight   Weighted timestep for this stage (e.g. for flux corrections)
subroutine ml_calcRHS(rhsType, rhsStruct, t, dtWeight)
   implicit none

   integer, intent(in) :: rhsType, rhsStruct
   real, intent(in) :: t
   real, intent(in) :: dtWeight

   return
end subroutine ml_calcRHS
