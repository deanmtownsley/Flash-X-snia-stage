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
!! @brief Spacetime_getLapse stub

!> @ingroup Spacetime
!!
!! @brief Obtain the lapse function at the provided location
!!
!! @details
!! @anchor Spacetime_getLapse_stub
!!
!! This procedure serves as an accessor to obtain the lapse function
!! at the specific location in the provided solution data.
!!
!! @param  alp       Scalar lapse function @f$ \alpha @f$
!! @param  solnData  Pointer to variables in UNK for the current tile
!! @param  loc       Location (i,j,k) in the current tile
subroutine Spacetime_getLapse(alp, solnData, loc)
#include "constants.h"

   implicit none

   real, intent(out) :: alp
   real, pointer :: solnData(:, :, :, :)
   integer, intent(in) :: loc(MDIM)

   alp = 0d0

   return
end subroutine Spacetime_getLapse
