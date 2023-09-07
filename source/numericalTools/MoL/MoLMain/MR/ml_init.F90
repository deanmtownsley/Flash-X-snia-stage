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
!! @brief ml_init implementation for MR

!> @ingroup MoLMR
!!
!! @brief Implements ml_init for MR
!!
!! @stubref{ml_init}
subroutine ml_init()
   use MoL_data, only: ml_nscratch
   use ml_mrData, only: FE, FI, FF, FAST_INITIAL, ml_slowMethod, ml_fastMethod, &
                        ml_nsubcycle, ml_nstages_fast, ml_nstages_slow

   use ml_interface, only: ml_initTableau

   use RuntimeParameters_interface, only: RuntimeParameters_get

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

   implicit none

   character(len=MAX_STRING_LENGTH) :: slowMethod_str, fastMethod_str

   integer :: i

   call RuntimeParameters_get("mr_slowMethod", slowMethod_str)
   call RuntimeParameters_get("mr_fastMethod", fastMethod_str)

   ml_slowMethod = trim(slowMethod_str)
   ml_fastMethod = trim(fastMethod_str)

   call RuntimeParameters_get("mr_nsubcycle", ml_nsubcycle)

   call ml_initTableau()

   ! Setup RHS indexing
   allocate (FE(ml_nstages_slow))
   allocate (FI(ml_nstages_slow))
   allocate (FF(ml_nstages_fast))

   FE = MOL_INVALID; FI = MOL_INVALID; FF = MOL_INVALID

   ! Uses MOL_RHS as first index
   do i = 1, ml_nstages_slow, 2
      FE(i) = MOL_RHS + i - 1
      FI(i) = MOL_RHS + i
   end do ! i

   FAST_INITIAL = FI(ml_nstages_slow - 1) + 1

   do i = 1, ml_nstages_fast
      FF(i) = FAST_INITIAL + i
   end do ! i

   ml_nscratch = ml_nstages_slow + ml_nstages_fast
end subroutine ml_init
