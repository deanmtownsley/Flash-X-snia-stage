!!****if* source/physics/sourceTerms/Flame/FlameSpeed/Constant/fl_flameSpeed
!!
!! NAME
!!
!!  fl_flameSpeed
!!
!! SYNOPSIS
!!
!!  call fl_flameSpeed(real, dimension(:,:,:,:),POINTER_INTENT_IN  :: solndata,
!!                     real, dimension(:,:,:)(out) :: flamespeed,
!!                     integer(in) :: blockid,
!!                     integer(in) :: nlayers)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!!
!! ARGUMENTS
!!
!!   solndata : solution data
!!
!!   flamespeed : flame speed
!!
!!   blockid : ID of block in current processor
!!
!!   nlayers : number of layers
!!
!!
!!
!!***


subroutine fl_flameSpeed( solnData, flamespeed, tileDesc, nlayers)

#include "Simulation.h"
#include "constants.h"
#include "FortranLangFeatures.fh"

  use fl_fsData, only : fl_fsConstFlameSpeed, fl_fsConstFlameWidth
  use Grid_tile,      ONLY : Grid_tile_t
  implicit none
  real, dimension(:,:,:,:),POINTER_INTENT_IN :: solnData
  real, dimension(:,:,:),intent(out) :: flamespeed
  type(Grid_tile_t), intent(in)     :: tileDesc
  integer, intent(in) :: nlayers

  flamespeed(:,:,:) = fl_fsConstFlameSpeed
#ifdef FSPD_VAR
  solndata(FSPD_VAR,:,:,:) = fl_fsConstFlameSpeed
#endif

end subroutine
