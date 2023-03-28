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
!! @brief ml_init implementation for ERK

!> @ingroup MoLERK
!!
!! @brief Implements ml_init for ERK
!!
!! @stubref{ml_init}
subroutine ml_init()
   use MoL_data, only: ml_nscratch
   use ml_erkData, only: ml_method, ml_stages, ml_K

   use ml_interface, only: ml_initTableau

   use RuntimeParameters_interface, only: RuntimeParameters_get

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

   implicit none

   character(len=MAX_STRING_LENGTH) :: erk_method_str

   integer :: i

   call RuntimeParameters_get("erk_method", erk_method_str)
   ml_method = trim(erk_method_str)

   call ml_initTableau

   ! Setup RHS indexing
   allocate (ml_K(ml_stages))

   ml_K = (/(MOL_RHS + i, i=0, ml_stages - 1)/)

   ml_nscratch = ml_stages - 1
end subroutine ml_init
