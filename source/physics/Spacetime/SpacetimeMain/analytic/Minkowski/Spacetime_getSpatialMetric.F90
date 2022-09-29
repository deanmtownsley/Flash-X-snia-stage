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
!! @brief Spacetime_getSpatialMetric implementation

!> @ingroup Minkowksi
!!
!! @brief Obtain the spatial metric at the provided location
!!
!! @details
!!
!! @stubref{Spacetime_getSpatialMetric}
subroutine Spacetime_getSpatialMetric(gxx, gxy, gxz, gyy, gyz, gzz, &
                                      tileDesc, solnData, loc)
   use Grid_tile, only: Grid_tile_t

#include "constants.h"

   implicit none

   real, intent(out) :: gxx, gxy, gxz, gyy, gyz, gzz
   type(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: solnData(:, :, :, :)
   integer, intent(in) :: loc(MDIM)

   gxx = 1d0
   gxy = 0d0
   gxz = 0d0
   gyy = 1d0
   gyz = 0d0
   gzz = 1d0
end subroutine Spacetime_getSpatialMetric
