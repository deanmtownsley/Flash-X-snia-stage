!!****if* source/physics/Hydro/HydroMain/Spark/Hydro
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
!!
!! NAME
!!
!!  Hydro
!!
!!
!! SYNOPSIS
!!
!!  Hydro(integer(IN) :: blockCount, <--deprecated
!!        integer(IN) :: blockList(blockCount) <--deprecated
!!        real(IN)    :: timeEndAdv,
!!        real(IN)    :: dt,
!!        real(IN)    :: dtOld,
!!        integer(IN) :: sweepOrder)
!!
!!
!! DESCRIPTION
!!
!!  Performs physics update in a directionally unsplit fashion.
!!
!!  The blockList and blockCount arguments tell this routine on
!!  which blocks and on how many to operate.  blockList is an
!!  integer array of size blockCount that contains the local
!!  block numbers of blocks on which to advance.
!!
!!  dt gives the timestep through which this update should advance,
!!  and timeEndAdv tells the time that this update will reach when
!!  it finishes.  dtOld gives the previously taken timestep.
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - end time
!!  dt         - timestep
!!  dtOld      - old timestep
!!  sweepOrder - dummy argument for the unsplit scheme, just a dummy
!!               variable to be consistent with a toplayer stub function
!!
!!***

!subroutine Hydro( nblk, blklst, &
!     timeEndAdv, dt,  dtOld,&
!     sweepOrder)

subroutine Hydro(timeEndAdv, dt, dtOld, sweepOrder)
  use Hydro_data, ONLY : hy_useHydro, hy_fluxCorrect, hy_gcMask, hy_lChyp, hy_C_hyp, hy_globalComm, hy_flx, hy_fly, hy_flz
  use hy_rk_interface, ONLY : hy_rk_eos, hy_rk_getFaceFlux, hy_rk_getGravAccel, hy_rk_updateSoln
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_conserveFluxes, &
      Grid_fillGuardCells, Grid_getTileIterator, &
      Grid_releaseTileIterator 
  use Eos_interface, ONLY : Eos_multiDim
  use IO_interface, ONLY : IO_setScalar
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t

  implicit none

#include "Simulation.h"
#include "constants.h"
  include "Flash_mpi.h"

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: blockDesc
  !integer              :: blockID
  !integer, intent(in) :: nblk !remove later
  !integer, intent(in) :: blklst(nblk) !remove later
  real,    intent(in) :: timeEndAdv, dt, dtOld
  integer, intent(IN) :: sweepOrder

  integer :: n, error
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(LOW:HIGH,MDIM) :: limits
  real :: hdt
  real, dimension(3) :: coeffs
  real, save :: onethird = 1./3.
  real, save :: twothird = 2./3.

  if (.NOT. hy_useHydro) return

  call Timers_start("Hydro")
  
  hdt = 0.5*dt

  ! Find the global maximum hyperbolic speed. hy_lChyp from Hydro_computeDt
#ifdef SPARK_GLM
  call MPI_AllReduce (hy_lChyp, hy_C_hyp, 1, &
       FLASH_REAL, MPI_MAX, hy_globalComm, error)
  call IO_setScalar("C_hyp", hy_lChyp)
#endif

  call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.false.,maskSize=NUNK_VARS,mask=hy_gcMask)

  !Mike Change to tile iterator here.  Check unsplit example
  ! Loop over blocks and compute Hydro update block-by-block
  !do n=1,nblk
  !   blockID = blklst(n)
  call Grid_getTileIterator(itor,LEAF,tiling=.TRUE.) !FALSE?
  do while(itor%isValid())
     call itor%currentTile(blockDesc)
     ! DivB will technically be lagged by 1 step, but we need ghost zones to
     ! compute the gradients. I ain't doing more communication for a diagnostic...
     call calcDivB(blockDesc)
     blkLimits(:,:)   = blockDesc%limits
     blkLimitsGC(:,:) = blockDesc%blkLimitsGC 
     blkLimits(:,1:NDIM) = blkLimits(:,1:NDIM)
     blkLimitsGC(:,1:NDIM) = blkLimitsGC(:,1:NDIM)
     call setLims(limits, NGUARD-1) !may remove this eventually
     call shockDetect(blockDesc,limits)
     ! Setup scratch storage of block data
     call Timers_start("scratch")
       ! U* = U0
     call saveState(blockDesc)
     call Timers_stop("scratch")
     !*********************
     ! Stage 1:
     !*********************

     ! calculate gravitational acceleration based on current value of GPOT_VAR
     ! This is stored in module-scope variable hy_grav
     call hy_rk_getGravAccel(blockDesc,limits)

     ! Perform reconstruction and flux calculation
     ! In Stage 1, compute low-side fluxes and update for NSTENCIL guardcells
     call setLims(limits, NSTENCIL)
     call Timers_start("getFaceFlux")
     call hy_rk_getFaceFlux(blockDesc, limits, hy_flx, hy_fly, hy_flz)
     call Timers_stop("getFaceFlux")
     if (hy_fluxCorrect) call addFluxes(0.5, .false.)
     ! Now update solution based on conservative fluxes
       ! U* = C1 * U0 + C2 * U* + C3 * dt*L(U*)
       ! U1 =  1 * U0           +  1 * dt*L(U0)
     coeffs = (/1.0, 0.0, 1.0/)
     call Timers_start("updateSoln")
     call hy_rk_updateSoln(blockDesc,dt,dtOld,limits,coeffs)
     call Timers_stop("updateSoln")
     ! Update EOS based on intermediate solution
    
     call Timers_start("eos")
     call hy_rk_eos(limits)
     call Timers_stop("eos")
     !*********************
     ! Stage 2:
     !*********************
     ! calculate gravitational acceleration based on current value of GPOT_VAR
     ! This is stored in module-scope variable hy_grav
     call hy_rk_getGravAccel(blockDesc,limits)
     
     ! Perform reconstruction and flux calculation
     ! In Stage 2, compute low-side fluxes and update for 0 guardcells
     call setLims(limits, 0)
     call Timers_start("getFaceFlux")
     call hy_rk_getFaceFlux(blockDesc, limits, hy_flx, hy_fly, hy_flz)
     call Timers_stop("getFaceFlux")
     if (hy_fluxCorrect) call addFluxes(0.5, .true.)
     ! Now update solution based on conservative fluxes
       ! U* =  C1 * U0 +  C2 * U* +  C3 * dt*L(U*)
       ! U2 = 1/2 * U0 + 1/2 * U1 + 1/2 * dt*L(U1)
     coeffs = (/0.5, 0.5, 0.5/)
     call Timers_start("updateSoln")
     call hy_rk_updateSoln(blockDesc,dt,dtOld,limits,coeffs)
     call Timers_stop("updateSoln")

     ! Update EOS based on intermediate solution
     call Timers_start("eos")
     call hy_rk_eos(limits)
     call Timers_stop("eos")
 
     ! Finally, store the output and free up the scratch array
     call Timers_start("scratch")
     call updateState(blockDesc)
     call Timers_stop("scratch")
     call itor%next()
  end do

  call Grid_releaseTileIterator(itor)
 
  print*, "Outside",hy_fluxCorrect
 
  if (hy_fluxCorrect) then
     ! Call the PARAMESH routine to correct fluxes at fine-coarse boundaries
     call Timers_start("flux correct")
     call Grid_conserveFluxes(ALLDIR,0)
     ! Loop over blocks and correct block-edge solution
     !do n=1,nblk
     !   blockID = blklst(n)
     call Grid_getTileIterator(itor,LEAF,tiling=.TRUE.) !FALSE?
     do while(itor%isValid())
        call itor%currentTile(blockDesc)
        blkLimits(:,:)   = blockDesc%limits
        blkLimitsGC(:,:) = blockDesc%blkLimitsGC
        !Mike this will require some work...
        call hy_rk_correctFluxes(blockDesc,dt)
        ! call Eos_multiDim(MODE_DENS_EI,blkLimits,blockID)
        call itor%next()
     end do
     call Timers_stop("flux correct")
     call Grid_releaseTileIterator(itor)
  end if

  ! Reset local maximum hyperbolic speed. This will be updated in Hydro_computeDt.
  hy_lChyp = TINY(1.0)

  call Timers_stop("Hydro")

contains

#include "Hydro_funcs.F90"

end subroutine Hydro
