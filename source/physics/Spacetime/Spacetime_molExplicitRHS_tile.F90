!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
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
!! @brief Spacetime_molExplicitRHS_tile stub

!> @ingroup physics_Spacetime
!!
!! @brief Adds Spacetime's contribution to the explicit RHS for a single tile
!!
!! @details
!! @anchor Spacetime_molExplicitRHS_tile_stub
!!
!! This procedure is responsible for adding its explicitly-integrated RHS
!! terms to the overall RHS for a single tile during the current integration stage
!!
!! @param  tileDesc   Descriptor for the current tile
!! @param  t          Time that the RHS is to be evaluated at
!! @param  dtWeight   Weighted timestep (e.g. for flux corrections)
subroutine Spacetime_molExplicitRHS_tile(tileDesc, t, dtWeight)
   use Grid_tile, only: Grid_tile_t

   implicit none

   type(Grid_tile_t), intent(in) :: tileDesc
   real, intent(in) :: t
   real, intent(in) :: dtWeight

   return
end subroutine Spacetime_molExplicitRHS_tile
