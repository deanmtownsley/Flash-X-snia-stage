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
!! @brief MoL_releaseFunctions implementation

!> @ingroup MoLMain
!!
!! @brief Implements MoL_releaseFunctions
!!
!! @stubref{MoL_releaseFunctions}
subroutine MoL_releaseFunctions()
   use ml_functions, only: ml_rhsE, ml_rhsI, ml_rhsF, &
                           ml_implicitUpdate, ml_postUpdate, ml_postUpdateFast, &
                           ml_rhsE_default, ml_rhsI_default, ml_rhsF_default, &
                           ml_implicitUpdate_default, ml_PostUpdate_default, &
                           ml_PostUpdateFast_default

   implicit none

   if (associated(ml_rhsE)) nullify (ml_rhsE)
   if (associated(ml_rhsI)) nullify (ml_rhsI)
   if (associated(ml_rhsF)) nullify (ml_rhsF)

   if (associated(ml_implicitUpdate)) nullify (ml_implicitUpdate)

   if (associated(ml_postUpdate)) nullify (ml_postUpdate)
   if (associated(ml_postUpdateFast)) nullify (ml_postUpdateFast)

   ml_rhsE => ml_rhsE_default
   ml_rhsI => ml_rhsI_default
   ml_rhsF => ml_rhsF_default

   ml_implicitUpdate => ml_implicitUpdate_default

   ml_postUpdate => ml_postUpdate_default
   ml_postUpdateFast => ml_postUpdateFast_default
end subroutine MoL_releaseFunctions
