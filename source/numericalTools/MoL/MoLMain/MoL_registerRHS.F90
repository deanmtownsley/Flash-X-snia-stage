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
!! @brief MoL_registerRHS implementation

!> @ingroup MoLMain
!!
!! @brief Implements MoL_registerRHS
!!
!! @stubref{MoL_registerRHS}
subroutine MoL_registerRHS(rhsType, rhsFunc)
   use MoL_functionTypes, only: MoL_rhs_t
   use ml_functions, only: ml_rhsE, ml_rhsI, ml_rhsF
   use ml_interface, only: ml_error

#include "MoL.h"

   implicit none

   integer, intent(in) :: rhsType
   procedure(MoL_rhs_t) :: rhsFunc

   select case (rhsType)
   case (MOL_RHS_EXPLICIT)
      if (associated(ml_rhsE)) nullify (ml_rhsE)
      ml_rhsE => rhsFunc

   case (MOL_RHS_IMPLICIT)
      if (associated(ml_rhsI)) nullify (ml_rhsI)
      ml_rhsI => rhsFunc

   case (MOL_RHS_FAST)
      if (associated(ml_rhsF)) nullify (ml_rhsF)
      ml_rhsF => rhsFunc

   case default
      call ml_error("Attempting to register unknown RHS function type")
   end select ! rhsType
end subroutine MoL_registerRHS
