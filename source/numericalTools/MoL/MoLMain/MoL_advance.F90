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
!! @brief MoL_advance implementation

!> @ingroup MoLMain
!!
!! @brief Implements MoL_advance
!!
!! @details
!! This implementation saves the intial state of the evolved variables at
!! time `t` and calls a method-specific `ml_advance`
!!
!! @stubref{MoL_advance}
subroutine MoL_advance(t, dt)
   use ml_interface, only: ml_advance
   use ml_memInterface, only: ml_memCopy

#include "MoL.h"

   implicit none

   real, intent(in) :: t, dt

   call ml_memCopy(MOL_INITIAL, MOL_EVOLVED)

   call ml_advance(t, dt)
end subroutine MoL_advance
