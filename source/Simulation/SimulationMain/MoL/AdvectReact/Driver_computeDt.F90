subroutine Driver_computeDt(nbegin, nstep, simTime, dtOld, dtNew)
    use Simulation_data

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile
    use Grid_iterator

#include "Simulation.h"
#include "constants.h"

    implicit none

    integer, intent(in)  :: nbegin, nstep
    real,    intent(in)  :: simTime    !! current simulation time
    real,    intent(in)  :: dtOld      !! last time step we used
    real,    intent(out) :: dtNew      !! the new timestep we get. to be returned.

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    real :: del(MDIM), dxMin

    dxMin = 1d99

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call tileDesc%deltas(del)

        dxMin = min(dxMin, del(IAXIS))

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)

    dtNew = sim_cfl * dxMin / sim_speed
end subroutine Driver_computeDt
