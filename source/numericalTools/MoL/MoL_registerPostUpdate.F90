!> @file source/numericalTools/MoL/MoL_registerPostUpdate.F90
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
!! @brief MoL_registerPostUpdate stub
!! @ingroup MoL

!> @brief Register a procedure responsible for performing post-update
!! @anchor MoL_registerPostUpdate_stub
!!
!! @param postUpdateType  post-update type identifier
!! @param postUpdateFunc  Procedure that will calculate the post update
!!
!! @pre `postUpdateType` is a valid MoL post-update identifier as defined in @ref MoL.h
!! @pre `postUpdateFunc` is a callable procedure of the current time: `call postUpdateFunc(t)`
!!
!! @sideffect MoL's procedure pointer for the specified post-update will be set to the
!!            provided procedure
!!
!! @returns None
!!
!! @details
!!    Valid post-update types include (defined in MoL.h):
!!       - `MOL_POST_UPDATE`      : Post-update (slow) per-stage
!!       - `MOL_POST_UPDATE_FAST` : Post-update (fast) per-stage
!!
!! @warning Will trigger Flash-X to abort if an invalid post-update type is specified
!!
!! @ingroup MoL
subroutine MoL_registerPostUpdate(postUpdateType, postUpdateFunc)
   use MoL_functionTypes, only: MoL_postUpdate_t

   implicit none

   integer, intent(in) :: postUpdateType
   procedure(MoL_postUpdate_t) :: postUpdateFunc

   return
end subroutine MoL_registerPostUpdate
