!> @file source/numericalTools/MoL/MoL_registerUpdate.F90
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
!! @brief MoL_registerUpdate stub
!! @ingroup MoL

!> @brief Register a procedure responsible for performing an update
!! @anchor MoL_registerUpdate_stub
!!
!! @param updateType  Update-type identifier
!! @param updateFunc  Procedure that will calculate the update
!!
!! @pre `updateType` is a valid MoL update identifier as defined in @ref MoL.h
!! @pre `updateFunc` is a callable procedure of the current time: `call updateFunc(t, dt)`
!!
!! @sideffect MoL's procedure pointer for the specified update will be set to the
!!            provided procedure
!!
!! @returns None
!!
!! @details
!!    Valid update types include (defined in MoL.h):
!!       - `MOL_IMPLICIT_UPDATE` : For equations and terms requiring implicit integration
!!
!! @note There is only one valid option for updateType (MOL_IMPLICIT_UPDATE),
!!       but this is left generic to accomodate for new update-types in the
!!       future, e.g. distinct implcit-updates for both slow- and fast-integration
!!       steps in the multi-rate integrator.
!!
!! @warning Will trigger Flash-X to abort if an invalid update-type is specified
!!
!! @ingroup MoL
subroutine MoL_registerUpdate(updateType, updateFunc)
   use MoL_functionTypes, only: MoL_update_t

   implicit none

   integer, intent(in) :: updateType
   procedure(MoL_update_t) :: updateFunc

   return
end subroutine MoL_registerUpdate
