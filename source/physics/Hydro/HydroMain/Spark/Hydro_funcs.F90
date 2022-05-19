!! *** source/physics/Hydro/HydroMain/Spark
!!
!! NAME
!! Hydro_funcs
!!
!! DESCRIPTION
!!
!! Holds functions frequently accessed by Hydro.F90 for Spark hydrodynamic solver.
!! 
!! *** 
!!Reorder directive used by FLASH with --index-reorder flag at setup
!!Reorder(4): hy_starState,solnData, U, hy_fl[xyz], hy_flxb[xyz]


subroutine addFluxes(weight,addFlux)
  !Store weighted fluxes, summed over RK stages, in temporary flux buffers.
  use Hydro_data, ONLY : hy_flx, hy_fly, hy_flz, hy_flxbx, hy_flxby, hy_flxbz 

  implicit none

#include "constants.h"
#include "Simulation.h"

  real, intent(IN) :: weight
  logical, intent(IN) :: addFlux
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: datasize
  integer :: lev, coarse, fine
  integer :: axis, geometry
  logical :: zeroFullRegister 
   
  lev = tileDesc%level
  blkLimits(:,:)   = tileDesc%limits

  if (addFlux) then
    hy_flxbx = hy_flxbx+weight*hy_flx(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
                                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) 
    if (NDIM > 1) &   
      hy_flxby = hy_flxby+weight*hy_fly(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,&
                                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
    if (NDIM > 2) &
      hy_flxbz = hy_flxbz+weight*hy_flz(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)
  else
    hy_flxbx = weight*hy_flx(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
                                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
    if (NDIM > 1) &   
      hy_flxby = weight*hy_fly(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,&
                                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
    if (NDIM > 2) &
      hy_flxbz = weight*hy_flz(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)
  endif
end subroutine addFluxes



!! Allocate variable size array holding local gravitational accelerations (depending on 
!! block size), fluxes, flux buffers, and save. 
!! If offloading to a device, send data and allocate on device only.
subroutine saveState(Uin,blkLimits,blkLimitsGC)
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Hydro_data, ONLY : hy_starState, hy_threadWithinBlock,hy_tmpState
  
  implicit none
#include "Simulation.h"
#include "constants.h"
#include "Spark.h"

  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real, pointer :: Uin(:,:,:,:)
  integer ::  v,i1,i2,i3
  ! call Timers_start("Allocations")

#ifdef OMP_OL
 !$omp target map(hy_starState, hy_tmpState)
 !$omp target teams distribute parallel do collapse(4) map(to:blkLimitsGC,Uin) &
 !$omp shared(Uin,hy_starState,blkLimitsGC,hy_tmpState) private(v,i1,i2,i3) default(none)
#else 
  !$omp parallel do collapse(4) if(hy_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared(Uin,hy_starState,blkLimits,blkLimitsGC,solState_tmp)
#endif
  do i3=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
    do i2=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
      do i1=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
        do v=1,NUNK_VARS
          hy_starState(v,i1,i2,i3) = Uin(v,i1,i2,i3)
          hy_tmpState(v,i1,i2,i3) = Uin(v,i1,i2,i3)
        enddo
      enddo
    enddo
  enddo

end subroutine saveState


subroutine updateState(Uin,blkLimits,blkLimitsGC)
  use Hydro_data, ONLY : hy_starState, hy_threadWithinBlock, hy_grav, hy_flx, hy_fly, hy_flz,&
                         hy_flxbx, hy_flxby, hy_flxbz, hy_fluxCorrect, hy_uplus, hy_uminus,&
                         hy_shk, hy_tposeBlk, hy_flux, hy_flat, hy_grv, hy_flat3d, hy_smalldens, hy_smallE, hy_smallpres, &
                         hy_smallX, hy_cvisc, hy_del, hy_tmpState
  implicit none
#include "Simulation.h"
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  real, pointer :: Uin(:,:,:,:)
#ifdef OMP_OL
!hy_uplus,hy_uminus,hy_shk,hy_tposeBlk,hy_flat,hy_grv,hy_flux
  !$omp target exit data map(DELETE:hy_flat3d,hy_tmpState,hy_flx,hy_fly,hy_flz)
  !$omp target exit data map(DELETE:hy_grav)
  !$omp target update from(hy_starState)
#endif
#ifndef FIXEDBLOCKSIZE
  if (hy_fluxCorrect) then
    deallocate(hy_flxbx);deallocate(hy_flxby);deallocate(hy_flxbz)
  endif
#endif

  !$omp parallel if(hy_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared(Uin,hy_starState,blkLimits,blkLimitsGC)
  !$omp workshare
#ifdef GPOT_VAR
  ! First reset GPOT_VAR.
  hy_starState(GPOT_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
       Uin(GPOT_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
#endif
  Uin(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
       hy_starState(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))

  !$omp end workshare
  !$omp end parallel
#ifdef OMP_OL
  !$omp target exit data map(DELETE:hy_starState)
#endif
  nullify(hy_starState)
  nullify(hy_tmpState)
end subroutine updateState

!! Set loop limits.  We include ngcell layers of guard zones
subroutine setLims(limits,ngcell,blkLimits)
  implicit none
#include "constants.h"
  integer, intent(IN):: ngcell
  integer, intent(IN),dimension(LOW:HIGH,MDIM) :: blkLimits
  integer, intent(OUT), dimension(LOW:HIGH,MDIM) :: limits
  integer :: ilim

  limits=blkLimits
  do ilim=1,NDIM
     limits(LOW ,ilim) = blkLimits(LOW ,ilim) - ngcell
     limits(HIGH,ilim) = blkLimits(HIGH,ilim) + ngcell
  end do
end subroutine setLims

!!****if* source/physics/Hydro/HydroMain/Spark/Hydro_funcs
!!
!! NAME
!!
!!  shockDetect
!!
!! SYNOPSIS
!!
!!  shockDetect( type(Grid_tile_t) :: tileDesc
!!               integer (IN)      :: limits )
!!
!! DESCRIPTION
!!
!!  This routine detects strongly compressive motions in simulation
!!  by calculating undivided pressure gradients and divergence of
!!  velocity fields. Two parameters beta and delta have been set
!!  to detect strong shocks. If such shocks exist then the unsplit
!!  scheme applies its robust flux differencings using limited slopes
!!  in data reconstruction step (see hy_rk_dataReconstruct.F90).
!!  Different shock strengths can also be detected by lowering/increasing
!!  beta and delta values.
!!
!! ARGUMENTS
!!
!!  tileDesc  - local block descriptor (respecting tiling syntax)
!!  limits     - region of the block in which to detect shocks
!!
!! REFERENCE
!!
!!  Balsara and Spicer, JCP, 149:270--292, 1999.
!!
!!***
subroutine shockDetect(tileDesc,limits)

  use Hydro_data,        only : hy_geometry, hy_tiny

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Spark.h"

  !! ---- Argument List ----------------------------------
  type(Grid_tile_t)   :: tileDesc
  integer, intent(IN) :: limits(LOW:HIGH,MDIM)
  !! -----------------------------------------------------

  integer :: i,j,k
  logical :: SW1, SW2

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

  real :: divv,gradPx,gradPy,gradPz
  real :: minP,minC,beta,delta
  real :: localCfl,cflMax
  real, dimension(:,:,:), allocatable :: Vc
  real, dimension(:,:,:,:), pointer   :: Uin
  
  !necessary for argument for %getDataPtr()
  nullify(Uin)

#ifndef SHOK_VAR
  return
#endif

  ! Two parameters that can be adjusted to detect shocks
  ! with different strengths:
  ! (a) The lower the values the weaker shocks detected
  !     (lower the values to detect more shock regions)
  ! (b) The larger the values the stronger shocks detected
  !     (increase the values to detect less shock regions)
  beta = 0.5 !0.5 !10. ! gradP
  delta= 0.1  !0.1 !2. ! divV
!!$  beta  = 0.1 !0.1
!!$  delta = 0.01


  blkLimits(:,:)   = tileDesc%limits
  blkLimitsGC(:,:) = tileDesc%blkLimitsGC

  call tileDesc%getDataPtr(Uin,CENTER)

  Uin(SHOK_VAR,:,:,:) = 0.
  
!!!*** keyword to indicate local allocation
  !! Allocate a temporary cell-centered array for sound speed
  allocate(Vc(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))

!!!*** keyword for GC loop nest
  !! Compute sound speed
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           Vc(i,j,k) = sqrt(Uin(GAMC_VAR,i,j,k)*Uin(PRES_VAR,i,j,k)&
                /max(Uin(DENS_VAR,i,j,k),hy_tiny))
        enddo
     enddo
  enddo

!!!*** keyword for non-standard loop nest
  do k=limits(LOW,KAXIS),limits(HIGH,KAXIS)
     do j=limits(LOW,JAXIS),limits(HIGH,JAXIS)
        do i=limits(LOW,IAXIS),limits(HIGH,IAXIS)

           ! initialize switch values
           SW1 = .false.
           SW2 = .false.

#if NDIM==1
           minP = minval(Uin(PRES_VAR,i-1:i+1,j,k))
           minC = minval(Vc(i-1:i+1,j,k))
#endif
#if NDIM==2
           minP = minval(Uin(PRES_VAR,i-1:i+1,j-1:j+1,k))
           minC = minval(Vc(i-1:i+1,j-1:j+1,k))
#endif
#if NDIM==3
           minP = minval(Uin(PRES_VAR,i-1:i+1,j-1:j+1,k-1:k+1))
           minC = minval(Vc(i-1:i+1,j-1:j+1,k-1:k+1))
#endif
           !! We do not need to include non-Cartesian geom factors here.
           !! Undivided divV
           divv =        Uin(VELX_VAR,i+1,j,  k  ) - Uin(VELX_VAR,i-1,j,  k  )
#if NDIM > 1
           divv = divv + Uin(VELY_VAR,i,  j+1,k  ) - Uin(VELY_VAR,i,  j-1,k  )
#if NDIM == 3
           divv = divv + Uin(VELZ_VAR,i,  j,  k+1) - Uin(VELZ_VAR,i,  j,  k-1)
#endif
#endif
           divv = 0.5*divv

           !! Undivided grad pres
           gradPx = 0.5*(Uin(PRES_VAR,i+1,j,  k  ) - Uin(PRES_VAR,i-1,j,  k  ))
           gradPy = 0.
           gradPz = 0.
#if NDIM > 1
           gradPy = 0.5*(Uin(PRES_VAR,i,  j+1,k  ) - Uin(PRES_VAR,i,  j-1,k  ))
#if NDIM == 3
           gradPz = 0.5*(Uin(PRES_VAR,i,  j,  k+1) - Uin(PRES_VAR,i,  j,  k-1))
#endif
#endif
           if ( abs(gradPx)+abs(gradPy)+abs(gradPz) .ge. beta*minP ) then
              SW1 = .true.
           endif
           if (-delta*minC .ge. divv) then
              SW2 = .true.
           endif
           if (SW1 .and. SW2) then
              ! Set SHOCK_VAR to 1.0 if a shock is detected.
              ! One use is for a local hybrid method in the Hydro unit which
              ! applies (a diffusive) HLL solver when SHOK_VAR = 1.
              Uin(SHOK_VAR,i,j,k) = 1.
           endif !endif (SW1 .and. SW2) then

        enddo !enddo i-loop
     enddo !enddo j-loop
  enddo !enddo k-loop

  ! Release block pointer
  call tileDesc%releaseDataPtr(Uin,CENTER)

  ! Deallocate sound speed array
  deallocate(Vc)

end subroutine shockDetect

!!Calculate divergence of the magnetic field.
subroutine calcDivB(tileDesc)
  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Spark.h"

  !! ---- Argument List ----------------------------------
  type(Grid_tile_t)     :: tileDesc
  !! -----------------------------------------------------

  integer :: i,j,k
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real, dimension(:,:,:,:), pointer   :: Uin
  real, dimension(MDIM) :: del
  real :: divB

  nullify(Uin)

#ifdef SPARK_GLM
  blkLimits(:,:)   = tileDesc%limits
  blkLimitsGC(:,:) = tileDesc%blkLimitsGC
  call tileDesc%getDataPtr(Uin,CENTER)
  call tileDesc%deltas(del)

!!!*** keyword for nonGC loop nest
  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           divB = 0.0
#if NDIM>1
           divB = (Uin(MAGX_VAR,i+1,j,k) - Uin(MAGX_VAR,i-1,j,k))&
                *0.5/del(IAXIS)
           divB = divB + (Uin(MAGY_VAR,i,j+1,k) - Uin(MAGY_VAR,i,j-1,k))&
                *0.5/del(JAXIS)
#if NDIM==3
           divB = divB + (Uin(MAGZ_VAR,i,j,k+1) - Uin(MAGZ_VAR,i,j,k-1))&
                *0.5/del(KAXIS)
#endif
#endif
           Uin(DIVB_VAR,i,j,k) = divB
        end do
     end do
  end do

  call tileDesc%releaseDataPtr(Uin,CENTER)
#endif
end subroutine calcDivB

