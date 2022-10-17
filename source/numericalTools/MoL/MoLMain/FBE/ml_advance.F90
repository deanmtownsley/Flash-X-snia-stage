!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! @brief ml_advance implementation for FBE

!> @ingroup MoLFBE
!!
!! @brief Implements ml_advance for FBE
!!
!! @stubref{ml_advance}
subroutine ml_advance(t, dt)
   use ml_functions, only: ml_postUpdate, ml_postUpdateFast, ml_implicitUpdate
   use ml_interface, only: ml_calcRHS
   use ml_memInterface, only: ml_memAddToVars

#include "MoL.h"

   implicit none

   real, intent(in) :: t, dt

   integer :: srcs(1)
   real :: facs(1)

   srcs = MOL_RHS
   facs = dt

   call ml_calcRHS(MOL_RHS_EXPLICIT, MOL_RHS, t)

   call ml_memAddToVars(MOL_EVOLVED, 1d0, 1, srcs, facs)

   call ml_postUpdate(t + dt)
   call ml_postUpdateFast(t + dt)

   call ml_implicitUpdate(t + dt, dt)
end subroutine ml_advance
