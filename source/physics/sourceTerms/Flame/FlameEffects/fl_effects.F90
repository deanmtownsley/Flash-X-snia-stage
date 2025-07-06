! Dean Townsley 2008
!
! this is a stub for flame effects that don't do anything locally

subroutine fl_effects( solnData, flamdot, dt, tileDesc)

  use Grid_tile, only: Grid_tile_t

  implicit none

  real, dimension(:,:,:,:),pointer,intent(in)  :: solnData
  real,dimension(:,:,:), intent(in)     :: flamdot
  real,intent(in)                       :: dt
  type(Grid_tile_t), intent(in)         :: tileDesc

  return
end subroutine
