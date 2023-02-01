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
!! @brief Spacetime_computeDt implementation

!> @ingroup Z4c
!!
!! @brief Implements Spacetime_computeDt for the Z4c solver
!!
!! @stubref{Spacetime_computeDt}
subroutine Spacetime_computeDt(tileDesc, solnData, dtMin, dtMinLoc)
   use Grid_tile, only: Grid_tile_t

   use z4c_data, only: z4c_cfl

#include "Simulation.h"
#include "constants.h"

   implicit none

   type(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: solnData(:, :, :, :)
   real, intent(inout) :: dtMin
   real, intent(inout) :: dtMinLoc(5)

   real, dimension(MDIM) :: del

   call tileDesc%deltas(del)

   ! The timestep required by the Z4c solver will always be
   ! limited by the speed of light (assumes c = 1), and is not
   ! tied to any specific location the tile
   dtMin = min(dtMin, z4c_cfl*minval(del(1:NDIM)))

end subroutine Spacetime_computeDt
