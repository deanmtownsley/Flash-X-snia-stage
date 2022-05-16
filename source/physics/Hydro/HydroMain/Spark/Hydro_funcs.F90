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
!!Reorder(4): hy_starState,solnData, U, hy_fl[xyz], hy_fluxBuf[XYZ]

subroutine select_RK_scheme(coeff_array,last_stage,limits_array,weights)
  !Select coefficients for update step, number of stages,
  !stencil sizes, and flux weights for each step.

  implicit none

#include "Simulation.h"

  real,intent(OUT), dimension(3,3) :: coeff_array
  integer, intent(OUT) :: last_stage
  integer, intent(OUT), dimension(3) :: limits_array
  real, intent(OUT), dimension(3) :: weights

  real, save :: onesixth = 1./6.
  real, save :: onethird = 1./3.
  real, save :: twothirds = 2./3.

#ifdef HY_RK3
  !RK3 quantities

  !Stage 1 coefficients
  ! U* = C1 * U0 + C2 * U* + C3 * dt*L(U*)
  ! U1 =  1 * U0           +  1 * dt*L(U0)

  !Stage 2 coefficients
  ! U* =  C1 * U0 +  C2 * U* +  C3 * dt*L(U*)
  ! U2 = 3/4 * U0 + 1/4 * U1 + 1/4 * dt*L(U1)

  !Stage 3 coefficients
  ! U* =  C1 * U0 +  C2 * U* +  C3 * dt*L(U*)
  ! U3 = 1/3 * U0 + 2/3 * U2 + 2/3 * dt*L(U2)
 
  !(remember FORTRAN is column major)
  coeff_array = reshape((/1.,0.75,onethird,0.,0.25,twothirds,1.,0.25,twothirds/),(/3,3/))

  last_stage = 3 

  !Array containing number of guard cells on each side for
  !the telescoping update.
  limits_array = (/2*NSTENCIL, NSTENCIL, 0/)

  !Weights that scale the fluxes as they're added into the buffers.
  !Here 'weights' is the same as 'coeff' used in Github pseudocode.
  weights = (/onesixth, onesixth, twothirds/)

#else
  !RK2 quantities

  ! Stage 1 coefficients
  ! U* = C1 * U0 + C2 * U* + C3 * dt*L(U*)
  ! U1 =  1 * U0           +  1 * dt*L(U0)

  ! Stage 2 coefficients
  ! Now update solution based on conservative fluxes
  ! U* =  C1 * U0 +  C2 * U* +  C3 * dt*L(U*)
  ! U2 = 1/2 * U0 + 1/2 * U1 + 1/2 * dt*L(U1)
  coeff_array = reshape((/1.,0.5,0.,0.,0.5,0.,1.,0.5,0./),(/3,3/))

  last_stage = 2

  limits_array = (/NSTENCIL, 0, 0/)

  weights = (/0.5,0.5,0./)
#endif

end subroutine select_RK_scheme

subroutine addFluxes(weight,addFlux)
  !Store weighted fluxes, summed over RK stages, in temporary flux buffers.
  use Hydro_data, ONLY : hy_flx, hy_fly, hy_flz, hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ 

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
    hy_fluxBufX = hy_fluxBufX+weight*hy_flx(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
                                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) 
    if (NDIM > 1) &   
      hy_fluxBufY = hy_fluxBufY+weight*hy_fly(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,&
                                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
    if (NDIM > 2) &
      hy_fluxBufZ = hy_fluxBufZ+weight*hy_flz(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)
  else
    hy_fluxBufX = weight*hy_flx(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
                                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
    if (NDIM > 1) &   
      hy_fluxBufY = weight*hy_fly(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,&
                                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
    if (NDIM > 2) &
      hy_fluxBufZ = weight*hy_flz(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                                           blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                                           blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1)
  endif
end subroutine addFluxes



!! Allocate variable size array holding local gravitational accelerations (depending on 
!! block size), fluxes, flux buffers, and save. 
!! If offloading to a device, send data and allocate on device only.
subroutine saveState(Uin,blkLimits,blkLimitsGC)
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Hydro_data, ONLY : hy_starState, hy_threadWithinBlock, hy_fluxCorrect, hy_grav, hy_flx, hy_fly, hy_flz,&
                         hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ, hydro_GPU_scratch, uPlusArray, uMinusArray,&
                         shck, snake, flux, flat, grv,hy_tiny,hy_hybridRiemann,hy_C_hyp,flat3d, &
                         hy_smalldens, hy_smallE, hy_smallpres, hy_smallX, hy_cvisc, del, solState_tmp,hy_geometry, &
                         hy_alphaGLM, scratch_allocated
  implicit none
#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS

  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real, pointer :: Uin(:,:,:,:)
  integer ::  v,i1,i2,i3
  ! call Timers_start("Allocations")

  ! Allocate needed space on GPU if it is not already there
  endif
#endif
!!$
!!$#ifndef FIXEDBLOCKSIZE
!!$  !Allocate size of flux buffers used for flux correction
!!$  if (hy_fluxCorrect) then
!!$    !allocate buffers here
!!$    if (.NOT. allocated(hy_fluxBufX)) then 
!!$      allocate(hy_fluxBufX(NFLUXES,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
!!$                                   blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
!!$                                   blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
!!$      hy_fluxBufX = 0.
!!$    endif
!!$
!!$    if (.NOT. allocated(hy_fluxBufY)) then 
!!$      allocate(hy_fluxBufY(NFLUXES,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
!!$                                   blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,&
!!$                                   blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
!!$      hy_fluxBufY = 0.
!!$    endif
!!$
!!$    if (.NOT. allocated(hy_fluxBufZ)) then 
!!$      allocate(hy_fluxBufZ(NFLUXES,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
!!$                                   blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
!!$                                   blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1))
!!$      hy_fluxBufZ = 0.
!!$    endif
!!$  endif
!!$  !Set up one block's worth of local gravity.  Allocation allows for compatibility with Paramesh4 and AMRex
!!$#endif
!!$  !Gravity 
!!$  if (.NOT. allocated(hy_grav)) then
!!$    allocate(hy_grav(MDIM,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
!!$      blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
!!$      blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
!!$  endif
!!$
!!$
!!$  if (.NOT. allocated(hy_starState)) then
!!$    allocate(hy_starState(NUNK_VARS,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
!!$      blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
!!$      blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
!!$  endif
!!$
!!$  if (.NOT. allocated(solState_tmp)) then
!!$    allocate(solState_tmp(NUNK_VARS,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
!!$      blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
!!$      blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
!!$  endif
!!$
!!$  ! update temp vars with solution data

#ifdef OMP_OL
 ! move data to GPU
 !$omp target enter data map(alloc:hy_starState,hy_flat3d,hy_tmpState,hy_flx,hy_fly,hy_flz,hy_grav)
 !$omp target update to(hy_tiny,hy_hybridRiemann,hy_C_hyp,hy_cvisc,hy_del,hy_smalldens, hy_smallE, hy_smallpres, hy_smallX,hy_geometry,hy_alphaGLM)
 ! distribute work throughout GPU
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


subroutine updateState(tileDesc)
  use Hydro_data, ONLY : hy_starState, hy_threadWithinBlock, hy_grav, hy_flx, hy_fly, hy_flz,&
                         hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ, hy_fluxCorrect, uPlusArray, uMinusArray,&
                         shck, snake, flux, flat, grv, flat3d, hy_smalldens, hy_smallE, hy_smallpres, &
                         hy_smallX, hy_cvisc, del, solState_tmp
  implicit none
#include "Simulation.h"
  type(Grid_tile_t),intent(IN) :: tileDesc
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real, pointer :: Uin(:,:,:,:)
#ifdef OMP_OL
!uPlusArray,uMinusArray,shck,snake,flat,grv,flux
  !$omp target exit data map(DELETE:flat3d,solState_tmp,hy_flx,hy_fly,hy_flz)
  !$omp target exit data map(DELETE:hy_grav)
  !$omp target update from(hy_starState)
#endif
  !First let's not forget to deallocate gravity
  deallocate(hy_grav)
  deallocate(hy_flx)
  deallocate(hy_fly)
  deallocate(hy_flz)
  deallocate(flat3d)
  deallocate(solState_tmp)

  ! deallocate(uPlusArray)
  ! deallocate(uMinusArray)
  ! deallocate(shck)
  ! deallocate(snake)
  ! deallocate(flat)
  ! deallocate(grv)
  ! deallocate(flux)

#ifndef FIXEDBLOCKSIZE
  if (hy_fluxCorrect) then
    deallocate(hy_fluxBufX);deallocate(hy_fluxBufY);deallocate(hy_fluxBufZ)
  endif
#endif
  nullify(Uin)
  call tileDesc%getDataPtr(Uin,CENTER)
 
  blkLimits(:,:)   = tileDesc%limits
  blkLimitsGC(:,:) = tileDesc%blkLimitsGC

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
  call tileDesc%releaseDataPtr(Uin,CENTER)
#ifdef OMP_OL
  !$omp target exit data map(DELETE:hy_starState)
#endif
  deallocate(hy_starState)
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

