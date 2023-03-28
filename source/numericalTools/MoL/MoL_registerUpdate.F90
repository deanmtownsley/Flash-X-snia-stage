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
!! @brief MoL_registerUpdate stub

!> @ingroup MoL
!!
!! @brief Register a procedure responsible for performing an update
!!
!! @details
!! @anchor MoL_registerUpdate_stub
!!
!! Valid update types include (defined in MoL.h):
!!    - `MOL_IMPLICIT_UPDATE` : For equations and terms requiring implicit integration
!!
!! @note There is only one valid option for updateType (MOL_IMPLICIT_UPDATE),
!!       but this is left generic to accomodate for new update-types in the
!!       future, e.g. distinct implcit-updates for both slow- and fast-integration
!!       steps in the multi-rate integrator.
!!
!! @param updateType  Update-type identifier
!! @param updateFunc  Procedure that will calculate the update
subroutine MoL_registerUpdate(updateType, updateFunc)
   use MoL_functionTypes, only: MoL_update_t

   implicit none

   integer, intent(in) :: updateType
   procedure(MoL_update_t) :: updateFunc

   return
end subroutine MoL_registerUpdate
