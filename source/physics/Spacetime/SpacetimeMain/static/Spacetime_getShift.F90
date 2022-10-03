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
!! @brief Spacetime_getShift implementation

!> @ingroup SpacetimeStatic
!!
!! @brief Obtain the shift vector at the provided location
!!
!! @details
!!
!! @stubref{Spacetime_getShift}
subroutine Spacetime_getShift(betax, betay, betaz, &
                              tileDesc, solnData, loc)
   use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "constants.h"

   implicit none

   real, intent(out) :: betax, betay, betaz
   type(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: solnData(:, :, :, :)
   integer, intent(in) :: loc(MDIM)

   betax = solnData(BETX_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
   betay = solnData(BETY_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
   betaz = solnData(BETZ_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
end subroutine Spacetime_getShift
