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
!! @brief Spacetime_getCurvature implementation

!> @ingroup SpacetimeMain
!!
!! @brief Obtain the extrinsic curvature at the provided location
!!
!! @details
!!
!! @stubref{Spacetime_getCurvature}
subroutine Spacetime_getCurvature(Kxx, Kxy, Kxz, Kyy, Kyz, Kzz, &
                                  tileDesc, solnData, loc)
   use Grid_tile, only: Grid_tile_t

#include "ADM.h"
#include "constants.h"

   implicit none

   real, intent(out) :: Kxx, Kxy, Kxz, Kyy, Kyz, Kzz
   type(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: solnData(:, :, :, :)
   integer, intent(in) :: loc(MDIM)

   Kxx = solnData(KXX_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
   Kxy = solnData(KXY_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
   Kxz = solnData(KXZ_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
   Kyy = solnData(KYY_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
   Kyz = solnData(KYZ_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
   Kzz = solnData(KZZ_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
end subroutine Spacetime_getCurvature
