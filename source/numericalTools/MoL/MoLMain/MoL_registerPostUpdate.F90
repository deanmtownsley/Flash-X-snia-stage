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
!! @brief MoL_registerPostUpdate implementation

!> @ingroup MoLMain
!!
!! @brief Implements MoL_registerPostUpdate
!!
!! @stubref{MoL_registerPostUpdate}
subroutine MoL_registerPostUpdate(postUpdateType, postUpdateFunc)
   use MoL_functionTypes, only: MoL_postUpdate_t
   use ml_functions, only: ml_postUpdate, ml_postUpdateFast
   use ml_interface, only: ml_error

#include "MoL.h"

   implicit none

   integer, intent(in) :: postUpdateType
   procedure(MoL_postUpdate_t) :: postUpdateFunc

   select case (postUpdateType)
   case (MOL_POST_UPDATE)
      if (associated(ml_postUpdate)) nullify (ml_postUpdate)
      ml_postUpdate => postUpdateFunc

   case (MOL_POST_UPDATE_FAST)
      if (associated(ml_postUpdateFast)) nullify (ml_postUpdateFast)
      ml_postUpdateFast => postUpdateFunc

   case default
      call ml_error("Attempting to register unknown post-update function type")
   end select ! postUpdateType
end subroutine MoL_registerPostUpdate
