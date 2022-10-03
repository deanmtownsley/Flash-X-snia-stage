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
!! @brief Spacetime_getLapse implementation

!> @ingroup SpacetimeStatic
!!
!! @brief Obtain the lapse function at the provided location
!!
!! @details
!!
!! @stubref{Spacetime_getLapse}
subroutine Spacetime_getLapse(alp, &
                              tileDesc, solnData, loc)
   use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "constants.h"

   implicit none

   real, intent(out) :: alp
   type(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: solnData(:, :, :, :)
   integer, intent(in) :: loc(MDIM)

   alp = solnData(ALP_VAR, loc(IAXIS), loc(JAXIS), loc(KAXIS))
end subroutine Spacetime_getLapse
