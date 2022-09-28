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
!! @brief Obtain all metric variables at the provided location
!!
!! @details
!! @anchor Spacetime_getMetric_stub
!!
!! This procedure serves as an accessor for all metric variables
!! at the specific location in the current grid tile.
!!
!! @param  metric    Contains all metric variables accesible via the
!!                   other accessors
!! @param  tileDesc  Descriptor for the current tile
!! @param  loc       Location (i,j,k) in the current tile
subroutine Spacetime_getMetric(metric, tileDesc, loc)
   use Spacetime_metric, only: Spacetime_metric_t
   use Grid_tile, only: Grid_tile_t

#include "constants.h"

   implicit none

   type(Spacetime_metric_t), intent(out) :: metric
   type(Grid_tile_t), intent(in) :: tileDesc
   integer, intent(in) :: loc(MDIM)

   ! Spacetime_metric_t has no data members at the stub-level

   return
end subroutine Spacetime_getMetric
