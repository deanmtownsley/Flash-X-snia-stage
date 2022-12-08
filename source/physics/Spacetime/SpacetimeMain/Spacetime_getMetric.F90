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
!! @brief Spacetime_getMetric implementation

!> @ingroup SpacetimeMain
!!
!! @brief Obtain the spatial metric at the provided location
!!
!! @details
!!
!! @stubref{Spacetime_getMetric}
subroutine Spacetime_getMetric(gxx, gxy, gxz, gyy, gyz, gzz, &
                               tileDesc, solnData, loc)
   use Grid_tile, only: Grid_tile_t

#include "ADM.h"
#include "constants.h"

   implicit none

   real, intent(out) :: gxx, gxy, gxz, gyy, gyz, gzz
   type(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: solnData(:, :, :, :)
   integer, intent(in) :: loc(MDIM)

   gxx = solnData(GXX_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
   gxy = solnData(GXY_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
   gxz = solnData(GXZ_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
   gyy = solnData(GYY_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
   gyz = solnData(GYZ_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
   gzz = solnData(GZZ_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
end subroutine Spacetime_getMetric
