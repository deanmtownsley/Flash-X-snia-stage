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
!! @brief Spacetime_getExtrinsicCurvature stub

!> @ingroup Spacetime
!!
!! @brief Obtain the extrinsic at the provided location
!!
!! @details
!! @anchor Spacetime_getExtrinsicCurvature_stub
!!
!! This procedure serves as an accessor to obtain the extrinsic curvature
!! at the specific location in the current grid tile.
!!
!! @param  Kxx,Kxy,Kxz,Kyy,Kyz,Kzz  Symmetric covariant components of the
!!                                  extrinsic @f$ K_{ij} @f$
!! @param  tileDesc                 Descriptor for the current tile
!! @param  solnData                 Pointer to variables in UNK for the current tile
!! @param  loc                      Location (i,j,k) in the current tile
subroutine Spacetime_getExtrinsicCurvature(Kxx, Kxy, Kxz, Kyy, Kyz, Kzz, &
                                           tileDesc, solnData, loc)
   use Grid_tile, only: Grid_tile_t

#include "constants.h"

   implicit none

   real, intent(out) :: Kxx, Kxy, Kxz, Kyy, Kyz, Kzz
   type(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: solnData(:, :, :, :)
   integer, intent(in) :: loc(MDIM)

   Kxx = 0d0
   Kxy = 0d0
   Kxz = 0d0
   Kyy = 0d0
   Kyz = 0d0
   Kzz = 0d0

   return
end subroutine Spacetime_getExtrinsicCurvature
