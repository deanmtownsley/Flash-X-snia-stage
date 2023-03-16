!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.

subroutine addFluxes(weight,addFlux)
  ! Now we handle storing the block face fluxes in PARAMESH arrays for flux correction
  ! Need to think about difference in size of flux arrays in PARAMESH and Spark
  ! Spark has more fluxes to account for mass scalars.  These do not need to exist in
  ! PARAMESH since they are just scaled mass fluxes...
  ! For now, Spark does not support flux correction. Conservation Smonservation.
  use Hydro_data, ONLY : hy_geometry, hy_fluxCorVars, &
       hy_flx, hy_fly, hy_flz
  use Grid_interface, ONLY : Grid_getBlkData, Grid_putFluxData

  implicit none

#include "constants.h"
#include "Simulation.h"

  real, intent(IN) :: weight
  logical, intent(IN) :: addFlux
  real, dimension(GRID_ILO_GC:GRID_IHI_GC,     &
       GRID_JLO_GC:GRID_JHI_GC,     &
       GRID_KLO_GC:GRID_KHI_GC) :: faceAreas

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: datasize

  blkLimits(:,:)   = blockDesc%limits
  blkLimitsGC(:,:) = blockDesc%blkLimitsGC
  
  datasize(1:MDIM) = blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

  if (hy_geometry /= CARTESIAN) then
     ! we are using consv_fluxes and need to divide by face areas
     call Grid_getBlkData(blockDesc,CELL_FACEAREA,ILO_FACE, EXTERIOR, &
          (/1,1,1/), faceAreas, datasize)
     call Grid_putFluxData(blockDesc%level,IAXIS,weight*hy_flx,datasize,hy_fluxCorVars,faceAreas,addFlux=addFLux)
     !print *, blockID, weight, faceAreas(7,1,1), hy_flx(1,7,1,1), faceAreas(27,1,1), hy_flx(1,27,1,1)
#if NDIM>1
     call Grid_getBlkData(blockDesc,CELL_FACEAREA,JLO_FACE, EXTERIOR, &
          (/1,1,1/), faceAreas, datasize)
     call Grid_putFluxData(blockDesc%level,JAXIS,weight*hy_fly,datasize,hy_fluxCorVars,faceAreas,addFlux=addFLux)
#if NDIM==3
     call Grid_getBlkData(blockDesc,CELL_FACEAREA,KLO_FACE, EXTERIOR, &
          (/1,1,1/), faceAreas, datasize)
     call Grid_putFluxData(blockDesc%level,KAXIS,weight*hy_flz,datasize,hy_fluxCorVars,faceAreas,addFlux=addFLux)
#endif
#endif
  else ! Cartesian geometry
     call Grid_putFluxData(blockDesc%level,IAXIS,weight*hy_flx,datasize,addFlux=addFLux)
     if (NDIM > 1) then
        call Grid_putFluxData(blockDesc%level,JAXIS,weight*hy_fly,datasize,addFlux=addFLux)
        if (NDIM > 2) then
           call Grid_putFluxData(blockDesc%level,KAXIS,weight*hy_flz,datasize,addFlux=addFLux)
        endif
     endif
  endif

end subroutine addFluxes


subroutine saveState(blockDesc)

  use Hydro_data, ONLY : hy_starState, hy_threadWithinBlock

  implicit none
#include "Simulation.h"

  type(Grid_tile_t),intent(IN)      :: blockDesc
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real, pointer :: solnData(:,:,:,:)
  integer :: xsize, ysize, zsize

  nullify(solnData)

  blkLimits(:,:)   = blockDesc%limits
  blkLimitsGC(:,:) = blockDesc%blkLimitsGC

  call blockDesc%getDataPtr(solnData,CENTER) 

  !old allocation strategy 
  !! Allocate storage for scratch array
  !xsize = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  !ysize = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  !zsize = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  !if (.NOT. allocated(hy_starState)) &
  !     allocate(hy_starState(NUNK_VARS,xsize,ysize,zsize))

  !new allocation allows for negative indices based on tiling strategy
  if (.NOT. allocated(hy_starState)) &
       allocate(hy_starState(NUNK_VARS,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))

  !$omp parallel if(hy_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared(solnData,hy_starState,blkLimits,blkLimitsGC)
  !$omp workshare
  hy_starState(:,:,:,:) = solnData(:,:,:,:)
  !$omp end workshare
  !$omp end parallel

  call blockDesc%releaseDataPtr(solnData,CENTER)

end subroutine saveState


subroutine updateState(blockDesc)
  use Hydro_data, ONLY : hy_starState, hy_threadWithinBlock
  implicit none
#include "Simulation.h"
  type(Grid_tile_t),intent(IN) :: blockDesc
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real, pointer :: solnData(:,:,:,:)
 
  nullify(solnData)
 
  blkLimits(:,:)   = blockDesc%limits
  blkLimitsGC(:,:) = blockDesc%blkLimitsGC

  call blockDesc%getDataPtr(solnData,CENTER)
  !$omp parallel if(hy_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared(solnData,hy_starState,blkLimits,blkLimitsGC)
  !$omp workshare
#ifdef GPOT_VAR
  ! First reset GPOT_VAR.
  hy_starState(GPOT_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
       solnData(GPOT_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
#endif
  solnData(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
       hy_starState(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
  !$omp end workshare
  !$omp end parallel
  call blockDesc%releaseDataPtr(solnData,CENTER)
  deallocate(hy_starState)
end subroutine updateState


subroutine setLims(limits,ngcell)
  implicit none
#include "constants.h"
  integer, intent(IN):: ngcell
  integer, intent(OUT), dimension(LOW:HIGH,MDIM) :: limits
  ! Set loop limits.  We include ngcell layers of guard zones
  limits = blkLimits
  limits(LOW ,IAXIS) = limits(LOW ,IAXIS) - ngcell
  limits(HIGH,IAXIS) = limits(HIGH,IAXIS) + ngcell
#if NDIM > 1
  limits(LOW ,JAXIS) = limits(LOW ,JAXIS) - ngcell
  limits(HIGH,JAXIS) = limits(HIGH,JAXIS) + ngcell
#if NDIM == 3
  limits(LOW ,KAXIS) = limits(LOW ,KAXIS) - ngcell
  limits(HIGH,KAXIS) = limits(HIGH,KAXIS) + ngcell
#endif /* NDIM == 3 */
#endif /* NDIM > 1 */
end subroutine setLims

!!****if* source/physics/Hydro/HydroMain/Spark/Hydro_funcs
!!
!! NAME
!!
!!  shockDetect
!!
!! SYNOPSIS
!!
!!  shockDetect( type(Grid_tile_t) :: blockDesc
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
!!  blockDesc  - local block descriptor (repsecting tiling syntax)
!!  limits     - region of the block in which to detect shocks
!!
!! REFERENCE
!!
!!  Balsara and Spicer, JCP, 149:270--292, 1999.
!!
!!***
subroutine shockDetect(blockDesc,limits)

  use Hydro_data,        only : hy_geometry, hy_tiny

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Spark.h"

  !! ---- Argument List ----------------------------------
  type(Grid_tile_t)   :: blockDesc
  integer, intent(IN) :: limits(LOW:HIGH,MDIM)
  !! -----------------------------------------------------

  integer :: i,j,k
  logical :: SW1, SW2

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

  real :: divv,gradPx,gradPy,gradPz
  real :: minP,minC,beta,delta
  real :: localCfl,cflMax
  real, dimension(:,:,:), allocatable :: Vc
  real, dimension(:,:,:,:), pointer   :: solnData
  
  !necessary for argument for %getDataPtr()
  nullify(solnData)

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


  blkLimits(:,:)   = blockDesc%limits
  blkLimitsGC(:,:) = blockDesc%blkLimitsGC
  call blockDesc%getDataPtr(solnData,CENTER)

  solnData(SHOK_VAR,:,:,:) = 0.

  !! Allocate a temporary cell-centered array for sound speed
  allocate(Vc(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))

  !! Compute sound speed
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           Vc(i,j,k) = sqrt(solnData(GAMC_VAR,i,j,k)*solnData(PRES_VAR,i,j,k)&
                /max(solnData(DENS_VAR,i,j,k),hy_tiny))
        enddo
     enddo
  enddo

  do k=limits(LOW,KAXIS),limits(HIGH,KAXIS)
     do j=limits(LOW,JAXIS),limits(HIGH,JAXIS)
        do i=limits(LOW,IAXIS),limits(HIGH,IAXIS)

           ! initialize switch values
           SW1 = .false.
           SW2 = .false.

#if NDIM==1
           minP = minval(solnData(PRES_VAR,i-1:i+1,j,k))
           minC = minval(Vc(i-1:i+1,j,k))
#endif
#if NDIM==2
           minP = minval(solnData(PRES_VAR,i-1:i+1,j-1:j+1,k))
           minC = minval(Vc(i-1:i+1,j-1:j+1,k))
#endif
#if NDIM==3
           minP = minval(solnData(PRES_VAR,i-1:i+1,j-1:j+1,k-1:k+1))
           minC = minval(Vc(i-1:i+1,j-1:j+1,k-1:k+1))
#endif
           !! We do not need to include non-Cartesian geom factors here.
           !! Undivided divV
           divv =        solnData(VELX_VAR,i+1,j,  k  ) - solnData(VELX_VAR,i-1,j,  k  )
#if NDIM > 1
           divv = divv + solnData(VELY_VAR,i,  j+1,k  ) - solnData(VELY_VAR,i,  j-1,k  )
#if NDIM == 3
           divv = divv + solnData(VELZ_VAR,i,  j,  k+1) - solnData(VELZ_VAR,i,  j,  k-1)
#endif
#endif
           divv = 0.5*divv

           !! Undivided grad pres
           gradPx = 0.5*(solnData(PRES_VAR,i+1,j,  k  ) - solnData(PRES_VAR,i-1,j,  k  ))
           gradPy = 0.
           gradPz = 0.
#if NDIM > 1
           gradPy = 0.5*(solnData(PRES_VAR,i,  j+1,k  ) - solnData(PRES_VAR,i,  j-1,k  ))
#if NDIM == 3
           gradPz = 0.5*(solnData(PRES_VAR,i,  j,  k+1) - solnData(PRES_VAR,i,  j,  k-1))
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
              solnData(SHOK_VAR,i,j,k) = 1.
           endif !endif (SW1 .and. SW2) then

        enddo !enddo i-loop
     enddo !enddo j-loop
  enddo !enddo k-loop

  ! Release block pointer
  call blockDesc%releaseDataPtr(solnData,CENTER)

  ! Deallocate sound speed array
  deallocate(Vc)

end subroutine shockDetect

subroutine calcDivB(blockDesc)
  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Spark.h"

  !! ---- Argument List ----------------------------------
  type(Grid_tile_t)     :: blockDesc
  !! -----------------------------------------------------

  integer :: i,j,k
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real, dimension(:,:,:,:), pointer   :: solnData
  real, dimension(MDIM) :: del
  real :: divB

  nullify(solnData)

#ifdef SPARK_GLM
  blkLimits(:,:)   = blockDesc%limits
  blkLimitsGC(:,:) = blockDesc%blkLimitsGC
  call blockDesc%getDataPtr(solnData,CENTER)
  call blockDesc%deltas(del)

  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           divB = 0.0
#if NDIM>1
           divB = (solnData(MAGX_VAR,i+1,j,k) - solnData(MAGX_VAR,i-1,j,k))&
                *0.5/del(IAXIS)
           divB = divB + (solnData(MAGY_VAR,i,j+1,k) - solnData(MAGY_VAR,i,j-1,k))&
                *0.5/del(JAXIS)
#if NDIM==3
           divB = divB + (solnData(MAGZ_VAR,i,j,k+1) - solnData(MAGZ_VAR,i,j,k-1))&
                *0.5/del(KAXIS)
#endif
#endif
           solnData(DIVB_VAR,i,j,k) = divB
        end do
     end do
  end do

  call blockDesc%releaseDataPtr(solnData,CENTER)
#endif
end subroutine calcDivB

!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_unitConvert
!!
!! NAME
!!
!!  hy_uhd_unitConvert
!!
!! SYNOPSIS
!!
!!  hy_uhd_unitConvert( integer (IN) :: blockID,
!!                      integer (IN) :: convertDir)
!!
!! DESCRIPTION
!!
!!  This routine converts physical unit of measure (CGS, SI, or NONE)
!!  for the unsplit MHD/Hydro units.
!!
!! ARGUMENTS
!!
!!  blockID    - local block ID
!!  convertDir - a direction of conversion (1 for forward, 0 for backward)
!!
!!***
subroutine unitConvert(blockDesc,convertDir)
  use Hydro_data,     ONLY : hy_bref

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Spark.h"

  !! Argument list -------------------------
  type(Grid_tile_t)   :: blockDesc
  integer, intent(IN) :: convertDir
  !! ---------------------------------------

  real, pointer, dimension(:,:,:,:) :: U

  nullify(U)

  call blockDesc%getDataPtr(U,CENTER)

  select case(convertDir)
  case(1)
#ifdef MAGX_VAR
     U(MAGX_VAR:MAGZ_VAR,:,:,:) = U(MAGX_VAR:MAGZ_VAR,:,:,:)/hy_bref
#endif

  case(0)
#ifdef MAGX_VAR
     U(MAGX_VAR:MAGZ_VAR,:,:,:) = U(MAGX_VAR:MAGZ_VAR,:,:,:)*hy_bref
#endif
  end select

  call blockDesc%releaseDataPtr(U,CENTER)
end subroutine unitConvert
