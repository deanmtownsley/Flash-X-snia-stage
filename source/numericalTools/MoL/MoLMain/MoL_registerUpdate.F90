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
!! @brief MoL_registerUpdate implementation

!> @ingroup MoLMain
!!
!! @brief Implements MoL_registerUpdate
!!
!! @stubref{MoL_registerUpdate}
subroutine MoL_registerUpdate(updateType, updateFunc)
   use MoL_functionTypes, only: MoL_update_t
   use ml_functions, only: ml_implicitUpdate
   use ml_interface, only: ml_error

#include "MoL.h"

   implicit none

   integer, intent(in) :: updateType
   procedure(MoL_update_t) :: updateFunc

   select case (updateType)
   case (MOL_IMPLICIT_UPDATE)
      if (associated(ml_implicitUpdate)) nullify (ml_implicitUpdate)
      ml_implicitUpdate => updateFunc

   case default
      call ml_error("Attempting to register unknown update function type")
   end select ! updateType
end subroutine MoL_registerUpdate
