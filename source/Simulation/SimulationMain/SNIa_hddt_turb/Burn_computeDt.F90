!! See source/physics/sourceTerms/Burn/Burn_computeDt.F90
!!  in flash distribution for main documention
!
! This version is the stub and should therefore not do anything.
! value of dt_burn is unmodified so it is not limiting.
!
! For type Ia simulations, we want to only enforce the hydro CFL,
! which is done elsewhere.  Although detonations can lead to
! burning which is unresolved (in time) with the hydrodynamic
! timestep, the interaction of PPM with a no-shock-burn condition
! leads to reasonable detonation speeds.  Flame burning, as
! implemented with the ADR flame tracking model is always on time
! scales longer than the hydrodynamic timestep for sub-sonic flame
! propagation speeds.
!
! Note that this makes certain presumptions about the burning
! (energy release source term) code's stability for timesteps which
! may be much larger than the burning timescale.
!
! Notes: Dean M. Townsley 2009

subroutine Burn_computeDt(blockID,                       & 
                           blkLimits,blkLimitsGC,        &
                           solnData,                     &
                           dt_burn, dt_minloc)


#include "constants.h"


  implicit none


  !! arguments
  integer, intent(IN)   :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real, pointer           :: solnData(:,:,:,:) 
  real, intent(INOUT)     :: dt_burn
  integer, intent(INOUT)  :: dt_minloc(5)

  return

end subroutine Burn_computeDt


