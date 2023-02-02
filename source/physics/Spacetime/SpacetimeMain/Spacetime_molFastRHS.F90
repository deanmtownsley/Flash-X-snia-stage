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
!! @brief Spacetime_molFastRHS implementation

!> @ingroup SpacetimeMain
!!
!! @brief Implements Spacetime_molFastRHS
!!
!! @stubref(Spacetime_molFastRHS)
subroutine Spacetime_molFastRHS(t, activeRHS, dtWeight)
   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator, Grid_fillGuardCells
   use Grid_iterator, only: Grid_iterator_t
   use Grid_tile, only: Grid_tile_t

   use Spacetime_interface, only: Spacetime_molFastRHS_tile

#include "constants.h"

   implicit none

   real, intent(in) :: t
   integer, intent(in) :: activeRHS
   real, intent(in) :: dtWeight

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   call Grid_fillGuardCells(CENTER, ALLDIR)

   call Grid_getTileIterator(itor, LEAF)

   TileLoop: do
      if (.not. itor%isValid()) exit TileLoop

      call itor%currentTile(tileDesc)

      call Spacetime_molFastRHS_tile(tileDesc, t, activeRHS, dtWeight)

      call itor%next()
   end do TileLoop

   call Grid_releaseTileIterator(itor)
end subroutine Spacetime_molFastRHS
