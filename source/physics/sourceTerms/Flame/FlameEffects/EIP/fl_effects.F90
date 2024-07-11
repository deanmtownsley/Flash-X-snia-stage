!!****if* source/physics/sourceTerms/Flame/FlameEffects/EIP/fl_effects
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
!!                  Grid_tile_t(in) :: tileDesc)
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!!
!! ARGUMENTS
!!
!!   solndata : 
!!
!!   flamdot : 
!!
!!   dt : 
!!
!!   tileDesc :
!!
!!
!!
!!***


subroutine fl_effects( solnData, flamdot, dt, tileDesc)

  use fl_effData, only : fl_effDeltae
  use Timers_interface, only : Timers_start, Timers_stop
  use Grid_interface, only : Grid_getBlkIndexLimits
  use Eos_interface, only : Eos_multiDim
  use Grid_tile, only: Grid_tile_t
  implicit none

#include "constants.h"
#include "Simulation.h"
#include "FortranLangFeatures.fh"
#include "Eos.h"
    
  real, dimension(:,:,:,:),POINTER_INTENT_IN  :: solnData
  real,dimension(:,:,:), intent(in)     :: flamdot
  real,intent(in)                       :: dt
  type(Grid_tile_t), intent(in)         :: tileDesc

  integer                    :: i, j, k
  real                       :: qdot
  integer,dimension(LOW:HIGH,MDIM) :: tileLimits

  tileLimits=tileDesc%limits

  ! update interior cells
  do k = tileLimits(LOW,KAXIS), tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS), tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS), tileLimits(HIGH,IAXIS)
              
           qdot = flamdot(i,j,k) * fl_effDeltae
           solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) + qdot*dt
#ifdef EINT_VAR
           solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + qdot*dt
#endif

        enddo
     enddo
  enddo
     
  ! We changed internal energy so need to update other quantities
  call Timers_start("eos")
  call Eos_multiDim(MODE_DENS_EI,tileLimits,solnData)
  !call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
  call Timers_stop("eos")
     
  return
end subroutine
