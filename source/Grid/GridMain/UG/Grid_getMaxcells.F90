subroutine Grid_getMaxcells(maxcells)
  use Grid_data, ONLY : gr_maxCells
  integer, intent(OUT) :: maxCells
  maxcells=gr_maxCells
end subroutine Grid_getMaxcells
