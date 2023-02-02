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
   use Driver_interface, only: Driver_getMype

   use Simulation_data, only: sim_dt

#include "constants.h"

   implicit none

   type(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: solnData(:, :, :, :)
   real, intent(inout) :: dtMin
   integer, intent(inout) :: dtMinLoc(5)

   integer :: myPE

   dtMin = sim_dt

   call Driver_getMype(MESH_COMM, myPE)

   dtMinLoc(1:3) = tileDesc%limits(LOW, :)
   dtMinLoc(4) = tileDesc%level
   dtMinLoc(5) = myPE

end subroutine Spacetime_computeDt
