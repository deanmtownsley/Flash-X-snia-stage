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
!! @brief ml_init implementation for IMEX

!> @ingroup MoLIMEX
!!
!! @brief Implements ml_init for IMEX
!!
!! @stubref{ml_init}
subroutine ml_init()
   use MoL_data, only: ml_nscratch
   use ml_imexData, only: ml_method, ml_nstages, FE, FI

   use ml_interface, only: ml_initTableau

   use RuntimeParameters_interface, only: RuntimeParameters_get

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

   implicit none

   character(len=MAX_STRING_LENGTH) :: imex_method_str

   integer :: i

   call RuntimeParameters_get("imex_method", imex_method_str)

   ml_method = trim(imex_method_str)

   call ml_initTableau()

   ! Setup RHS indexing
   allocate (FE(ml_nstages))
   allocate (FI(ml_nstages))

   FE = MOL_INVALID; FI = MOL_INVALID

   ! Uses MOL_RHS as first index
   do i = 0, ml_nstages - 1
      FE(i + 1) = MOL_RHS + 2*i
      FI(i + 1) = MOL_RHS + 2*i + 1
   end do ! i

   ! -1 for reusing MOL_RHS
   ml_nscratch = 2*ml_nstages - 1
end subroutine ml_init
