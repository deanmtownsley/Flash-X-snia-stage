#include "Simulation.h"
subroutine Grid_getMaxcells(maxcells)
  integer, intent(OUT) :: maxCells
  maxcells=max(NXB,NYB,NZB)
end subroutine Grid_getMaxcells
