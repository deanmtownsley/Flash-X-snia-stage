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
!! @brief Spacetime_getMetric stub

!> @ingroup Spacetime
!!
!! @brief Obtain the spatial metric at the provided location
!!
!! @details
!! @anchor Spacetime_getMetric_stub
!!
!! This procedure serves as an accessor to obtain the spatial metric
!! at the specific location in the provided solution data.
!!
!! @param  gxx,gxy,gxz,gyy,gyz,gzz  Symmetric covariant components of the
!!                                  spatial metric @f$ g_{ij} @f$
!! @param  tileDesc                 Current tile descriptor
!! @param  solnData                 Pointer to variables in UNK for the current tile
!! @param  loc                      Location (i,j,k) in the current tile
subroutine Spacetime_getMetric(gxx, gxy, gxz, gyy, gyz, gzz, &
                               tileDesc, solnData, loc)
   use Grid_tile, only: Grid_tile_t

#include "constants.h"

   implicit none

   real, intent(out) :: gxx, gxy, gxz, gyy, gyz, gzz
   type(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: solnData(:, :, :, :)
   integer, intent(in) :: loc(MDIM)

   gxx = 0d0
   gxy = 0d0
   gxz = 0d0
   gyy = 0d0
   gyz = 0d0
   gzz = 0d0

   return
end subroutine Spacetime_getMetric
