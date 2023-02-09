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
!! @brief Spacetime_molPostTimeStep implementation

!> @ingroup Z4c
!!
!! @brief Implements Spacetime_molPostTimeStep for the Z4c solver
!!
!! @stubref{Spacetime_molPostTimeStep}
subroutine Spacetime_molPostTimeStep(t)
   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator, Grid_fillGuardCells
   use Grid_iterator, only: Grid_iterator_t
   use Grid_tile, only: Grid_tile_t

   use z4c_interface, only: z4c_calculateConstraintViolation

#include "Z4c.h"
#include "constants.h"

   implicit none

   real, intent(in) :: t

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   integer :: i, j, k, lim(LOW:HIGH, MDIM)
   real, pointer :: vars(:, :, :, :)
   real :: del(MDIM)

   nullify (vars)

   call Grid_fillGuardCells(CENTER, ALLDIR)

   call Grid_getTileIterator(itor, LEAF)

   TileLoop: do
      if (.not. itor%isValid()) exit TileLoop

      call itor%currentTile(tileDesc)

      lim = tileDesc%limits
      lim(LOW, IAXIS) = lim(LOW, IAXIS) + 4
      lim(HIGH, IAXIS) = lim(HIGH, IAXIS) - 4

      call tileDesc%deltas(del)

      call tileDesc%getDataPtr(vars, CENTER)

      call z4c_calculateConstraintViolation(vars, lim, del)

      call tileDesc%releaseDataPtr(vars, CENTER)

      call itor%next()
   end do TileLoop

   call Grid_releaseTileIterator(itor)
end subroutine Spacetime_molPostTimeStep
