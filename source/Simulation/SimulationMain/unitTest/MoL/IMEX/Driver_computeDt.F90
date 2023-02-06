subroutine Driver_computeDt(nbegin, nstep, simTime, dtOld, dtNew)
   use Simulation_data, only: sim_alpha

   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator, &
                             Grid_getCellCoords
   use Grid_iterator, only: Grid_iterator_t
   use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "constants.h"

   implicit none

   integer, intent(in)  :: nbegin, nstep
   real, intent(in)  :: simTime    !! current simulation time
   real, intent(in)  :: dtOld      !! last time step we used
   real, intent(out) :: dtNew      !! the new timestep we get. to be returned.

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   real, dimension(MDIM) :: del
   real :: dtMin

   dtMin = huge(1d0)

   call Grid_getTileIterator(itor, LEAF)

   TileLoop: do
      if (.not. itor%isValid()) exit TileLoop

      call itor%currentTile(tileDesc)

      call tileDesc%deltas(del)

      dtMin = min(dtMin, 0.5d0*del(IAXIS)/sim_alpha)

      call itor%next()
   end do TileLoop

   call Grid_releaseTileIterator(itor)

   dtNew = dtMin
end subroutine Driver_computeDt
