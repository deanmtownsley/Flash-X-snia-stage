!!****if* source/physics/sourceTerms/Flame/FlameEffects/BurnParametric/fl_effects
!!
!! NAME
!!
!!  fl_effects
!!
!! SYNOPSIS
!!
!!  call fl_effects(real, dimension(:,:,:,:),POINTER_INTENT_IN  :: solndata,
!!                  real,dimension(:,:,:)(in) :: flamdot,
!!                  real(in) :: dt,
!!                  integer(in) :: blockid)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!!
!! save the time derivative representing the actual change in the
!! flam scalar for use by the burn source term
!! TODO should change FLDT to a scratch variable once things have stabilized
!!
!! ARGUMENTS
!!
!!   solndata : 
!!
!!   flamdot : 
!!
!!   dt : 
!!
!!   blockid : ID of block in current processor
!!
!!
!!
!!***



#include "Simulation.h"
#include "constants.h"
#include "FortranLangFeatures.fh"
subroutine fl_effects( solnData, flamdot, dt, tileDesc)

  use Grid_tile, only: Grid_tile_t

  implicit none

  real, dimension(:,:,:,:),pointer,intent(in)      :: solnData
  real,dimension(:,:,:), intent(in), allocatable   :: flamdot
  real,intent(in)                                  :: dt
  type(Grid_tile_t), intent(in)                    :: tileDesc

  ! only need interior cells
  solnData(FLDT_VAR,tileDesc%limits(LOW,IAXIS):tileDesc%limits(HIGH,IAXIS),   &
                    tileDesc%limits(LOW,JAXIS):tileDesc%limits(HIGH,JAXIS),   &
                    tileDesc%limits(LOW,KAXIS):tileDesc%limits(HIGH,KAXIS)) = &
                         flamdot(tileDesc%limits(LOW,IAXIS):tileDesc%limits(HIGH,IAXIS),   &
                                 tileDesc%limits(LOW,JAXIS):tileDesc%limits(HIGH,JAXIS),   &
                                 tileDesc%limits(LOW,KAXIS):tileDesc%limits(HIGH,KAXIS))

  return
end subroutine
