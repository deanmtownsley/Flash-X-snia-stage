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
!! @brief Spacetime_computeDt stub

!> @ingroup Spacetime
!!
!! @brief Calculate the minimum required timestep
!!
!! @details
!! @anchor Spacetime_computeDt_stub
!!
!! This procedure will calculate the minimum necessary timestep
!! required by the Spacetime unit
!!
!! @param  tileDesc  Information and access to the solution data for
!!                   the current tile
!! @param  dtMin     On output, the minimum required timstep that
!!                   is required by the Spacetime unit if it is
!!                   smaller than the input value
!! @param  dtMinLoc  The location in the grid responsible for the
!!                   minumum required timestep, indexed as
!!                   `(i,j,k,block,proc)`
subroutine Spacetime_computeDt(tileDesc, dtMin, dtMinLoc)
   use Grid_tile, only: Grid_tile_t

   implicit none

   type(Grid_tile_t), intent(in) :: tileDesc
   real, intent(inout) :: dtMin
   real, intent(inout) :: dtMinLoc(5)

   return
end subroutine Spacetime_computeDt
