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
!!
!! NAME
!!
!!  Hydro
!!
!!   For more details see the documentation of the NULL implementation
!!
!!***
!!Reorder(4): hy_fl[xyz]

subroutine Hydro(timeEndAdv, dt, dtOld, sweepOrder)

  use Hydro_data, ONLY :   hya_starState, &      
       hya_flx, hya_fly, hya_flz, hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ,&
       hy_farea,hy_cvol,hy_xCenter,hy_xLeft,hy_xRight,hy_yCenter,hy_zCenter
  
  use Hydro_data, ONLY : hy_fluxCorrect, hy_fluxCorrectPerLevel,hy_gcMask,&
       hy_lChyp, hy_C_hyp, hy_geometry,hy_del, hy_tiny, hy_hybridRiemann

  use hy_rk_interface, ONLY : hy_rk_getFaceFlux, hy_rk_getGraveAccel, hy_rk_updateSoln, hy_rk_correctFluxes
  use hy_data, ONLY : hy_useHydro
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_communicateFluxes, &
       Grid_fillGuardCells, Grid_getMaxRefinement, &
       Grid_correctFluxData_xtra, Grid_putFluxData_block, &
       Grid_getFluxCorrData_block,Grid_renormAbundance
  use Grid_interface, ONLY : Grid_getCellCoords, Grid_getCellFaceAreas, &
       Grid_getCellVolumes
  
  use Eos_interface, ONLY : Eos_wrapped
  use IO_interface, ONLY : IO_setScalar
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator
  use Grid_tile,         ONLY : Grid_tile_t
  use Grid_iterator,     ONLY : Grid_iterator_t
  
  implicit none
  
#include "Simulation.h"
#include "constants.h"
#include "Eos.h"
#include "Spark.h"
  include "Flashx_mpi.h"
  
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS
  type(Grid_iterator_t) :: itor
  real,              pointer    :: Uin(:,:,:,:)
  type(Grid_tile_t)     :: tileDesc
  integer :: level
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC,grownLimits
  real,dimension(MDIM) :: deltas

  real,    intent(in) :: timeEndAdv, dt, dtOld
  integer, intent(IN) :: sweepOrder
  integer :: n, error, maxLev=-1
  integer :: xLo,xHi,yLo,yHi,zLo,zHi
  integer, dimension(LOW:HIGH,MDIM) :: limits
  real :: hdt
  integer :: stage, last_stage
  real, dimension(3) :: coeffs, weights
  integer, dimension(3) :: limits_array
  real, dimension(3,3) :: coeff_array
    logical, dimension (3) :: addFlux_array
  integer :: i,j,k,v,maxcells
  
    real, pointer,dimension(:,:,:,:) :: hy_starState,hy_flx,hy_fly,hy_flz
  
  integer, dimension(MDIM) :: lo, hi, loGC, hiGC
  integer :: xLoGC,yLoGC,zLoGC,xHiGC,yHiGC,zHiGC
  integer :: pLo,pHi !low and high indices for pencil arrays
  integer :: lev
  
  if (.NOT. hy_useHydro) return
  
  call Timers_start("Hydro")
  
  call check_if_on_GPU()

  
  hdt = 0.5*dt
  
  !Array indicating whether to add flux into flux buffers (True)
  !or overwrite them (False).
  !This array will work for RK2 & RK3
  addFlux_array = (/.false.,.true.,.true./)
  
  !set up quantities specific to RK scheme (lives in Hydro_funcs)

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
  !Weights that scale the fluxes as they are added into the buffers.
  !Here weights is
  the same as coeff used in Github pseudocode.
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
  
  ! Find the global maximum hyperbolic speed. hy_lChyp from Hydro_computeDt
#ifdef SPARK_GLM
  call MPI_AllReduce (hy_lChyp, hy_C_hyp, 1, &
       FLASH_REAL, MPI_MAX, hy_globalComm, error)
  call IO_setScalar("C_hyp", hy_lChyp)
#endif
  
#ifdef FLASH_GRID_UG
  hy_fluxCorrect = .false.
  maxLev = 1
#else  
  ! mode=1 means lrefine_max, which does not change during sim.
  call Grid_getMaxRefinement(maxLev, mode=1)
#endif
  
  call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.false.,maskSize=NUNK_VARS,mask=hy_gcMask)
  !-------------------------------------------------------------------!
  !***NO Flux correction    or   Flux correction but NOT per level****!
  !-------------------------------------------------------------------!
  if ((.NOT.hy_fluxCorrect).OR.((hy_fluxCorrect).AND.(.NOT.hy_fluxCorrectPerLevel))) then
     ! Loop over blocks and compute Hydro update block-by-block
     nullify(Uin)
     call Grid_getTileIterator(itor, LEAF, tiling=.false.)
     do while(itor%isValid())
        call itor%currentTile(tileDesc)
        blkLimits(:,:)=tileDesc%limits
        blkLimitsGC(:,:)=tileDesc%blkLimitsGC
        grownLimits(:,:)=tileDesc%grownLimits
        call tileDesc%deltas(deltas)
        level=tileDesc%level
        call tileDesc%getDataPtr(Uin, CENTER)  
     ! DivB will technically be lagged by 1 step, but we need ghost zones to
     ! compute the gradients. I ain't doing more communication for a diagnostic...
        lo(:) = blkLimits(LOW,:)
        hi(:) = blkLimits(HIGH,:)
        loGC(:) = blkLimitsGC(LOW,:)
        hiGC(:) = blkLimitsGC(HIGH,:)
        !convenience indices
        xLoGC = loGC(IAXIS); xHiGC = hiGC(IAXIS)
        yLoGC = loGC(JAXIS); yHiGC = hiGC(JAXIS)
        zLoGC = loGC(KAXIS); zHiGC = hiGC(KAXIS)
        
        call allocate_scr(blkLimits,blkLimitsGC)
        if (hy_geometry /= CARTESIAN) then
           allocate(hy_farea(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC))
           call Grid_getCellFaceAreas(IAXIS,level,loGC,hiGC,hy_farea)
           allocate(hy_cvol(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC))
           call Grid_getCellVolumes(level,loGC,hiGC,hy_cvol)
           allocate(hy_xCenter(xLoGC:xHiGC))
           call Grid_getCellCoords(IAXIS, CENTER, level, loGC, hiGC, hy_xCenter)
           allocate(hy_xLeft(xLoGC:xHiGC))
           call Grid_getCellCoords(IAXIS, LEFT_EDGE, level, loGC, hiGC, hy_xLeft)
           allocate(hy_xRight(xLoGC:xHiGC))
           call Grid_getCellCoords(IAXIS, RIGHT_EDGE, level, loGC, hiGC, hy_xRight)
           allocate(hy_yCenter(yLoGC:yHiGC))
           call Grid_getCellCoords(JAXIS, CENTER, level, loGC, hiGC, hy_yCenter)
           allocate(hy_zCenter(zLoGC:zHiGC))
           call Grid_getCellCoords(KAXIS, CENTER, level, loGC, hiGC, hy_zCenter)
        endif
        hy_del=deltas
        call setLims(NGUARD-1,blkLimits,limits)
     
        !---JaredC--- Action 2 (no dependence within interation) 
        nullify(Uin)
        call tileDesc%getDataPtr(Uin,CENTER)
     
        call calcDivB(Uin,hy_del,blkLimits) !---JaredC--- This could be done on GPU or CPU and asynchronous to the rest of the sections
        !---JaredC--- This could be pulled out of the loop and given its own loop.
     
        call shockDetect(Uin,limits,blkLimitsGC)

        ! Setup scratch storage of block data (grav included)
        !---JaredC--- Send data to device and/or allocate arrays
        call Timers_start("Offloaded Section")
        call Timers_start("scratch")
        ! U* = U0
        !  The data is sent to the GPU at this point.
        call saveState(Uin,blkLimits,blkLimitsGC)
        
        
        call Timers_stop("scratch")
        
        !---JaredC--- Action 3 
        !Begin loop over stages
        do stage=1,last_stage
           ! calculate gravitational acceleration based on current value of GPOT_VAR
           ! This is stored in module-scope variable hy_grav
           call  hy_rk_getGraveAccel (hy_del,limits,blkLimitsGC)
           
           !Set needed number of guard cells to update based on
           !current stage for telescoping upda
           call setLims(limits_array(stage), blkLimits, limits)
           
           call Timers_start("getFaceFlux")
           
           call hy_rk_getFaceFlux(blklimits,blkLimitsGC, limits)
           

           call Timers_stop("getFaceFlux")
           !------- Add this in -----
           if (hy_fluxCorrect) call addFluxes(level,blkLimits,weights(stage),addFlux_array(stage)) 
           
           ! Now update solution based on conservative fluxes
           ! See select_RK_scheme() for explicit outline of C1, C2, & C3
           ! U* = C1 * U0 + C2 * U* + C3 * dt*L(U*)
           
           !Set proper coefficients for the given stage
           coeffs = coeff_array(stage,:)
           
           call Timers_start("updateSoln")
           call hy_rk_updateSoln(Uin,blkLimits,blklimitsGC,level,hy_del, dt, dtOld, limits, coeffs)
           call Timers_stop("updateSoln")
           
#if NSPECIES>0
           !Properly normalize species after the update
           call Driver_abort("Grid_renormAbundance not implemented in SPARK with GPU offloading yet")
#endif 
           
           ! Update EOS based on intermediate solution
           call Timers_start("eos")
           hy_starState(1:NUNK_VARS,&
                blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))=>hya_starState
           call Eos_wrapped(MODE_DENS_EI,limits,hy_starState)
           call Timers_stop("eos")
           nullify(hy_starState)
           !Finally update the state
           if (stage == last_stage) then
              call Timers_start("scratch")
              call updateState(Uin,blkLimits,blkLimitsGC)
              call Timers_stop("scratch")
           endif
           
        enddo!stage loop
        if (hy_geometry /= CARTESIAN) then
           deallocate(hy_xCenter)
           deallocate(hy_xLeft)
           deallocate(hy_xRight)
           deallocate(hy_farea)
           deallocate(hy_cvol)
           deallocate(hy_yCenter)
           deallocate(hy_zCenter)
        end if
        !Store flux buffer in semipermanent flux storage (SPFS) 
        if (hy_fluxCorrect) call Grid_putFluxData_block(tileDesc,&
             hy_fluxBufX,hy_fluxBufY,hy_fluxBufZ,blkLimits(LOW,:)) 
        call Timers_stop("Offloaded Section")
        call deallocate_scr

        call tileDesc%releaseDataPtr(Uin,CENTER)
        call itor%next()
     end do !!block loop
     call Grid_releaseTileIterator(itor)
        
     
     if (hy_fluxCorrect) then
        !Communicate the fine fluxes
        call Grid_communicateFluxes(ALLDIR,UNSPEC_LEVEL)
        
        !        call Timers_start("flux correct")
        !        ! Loop over blocks and correct block-edge solution
        
        nullify(Uin)
        call Grid_getTileIterator(itor, LEAF, tiling=.false.)
        do while(itor%isValid())
           call itor%currentTile(tileDesc)
           blkLimits(:,:)=tileDesc%limits
           blkLimitsGC(:,:)=tileDesc%blkLimitsGC
           grownLimits(:,:)=tileDesc%grownLimits
           call tileDesc%deltas(deltas)
           level=tileDesc%level
           call tileDesc%getDataPtr(Uin, CENTER) 
           hy_del=deltas
        
           
           !           !Get 'Flux hy_del' on coarse side of fine coarse boundaries; 
           !           !all other values are 0.
           call allocate_fxscr(blkLimits,blkLimitsGC)
           
           call Grid_getFluxCorrData_block(tileDesc,hy_fluxBufX,hy_fluxBufY,hy_fluxBufZ,&
                blkLimits(LOW,:))
           
           if (hy_geometry /= CARTESIAN) then
              allocate(hy_farea(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC))
              call Grid_getCellFaceAreas(IAXIS,level,loGC,hiGC,hy_farea)
              allocate(hy_cvol(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC))
              call Grid_getCellVolumes(level,loGC,hiGC,hy_cvol)
              allocate(hy_xCenter(xLoGC:xHiGC))
              call Grid_getCellCoords(IAXIS, CENTER, level, loGC, hiGC, hy_xCenter)
              allocate(hy_xLeft(xLoGC:xHiGC))
              call Grid_getCellCoords(IAXIS, LEFT_EDGE, level, loGC, hiGC, hy_xLeft)
              allocate(hy_xRight(xLoGC:xHiGC))
              call Grid_getCellCoords(IAXIS, RIGHT_EDGE, level, loGC, hiGC, hy_xRight)
           endif
        
           
           call hy_rk_correctFluxes(Uin,blkLimits,blklimitsGC,level,hy_del, dt)
           if (hy_geometry /= CARTESIAN) then
              deallocate(hy_xCenter)
              deallocate(hy_xLeft)
              deallocate(hy_xRight)
              deallocate(hy_farea)
              deallocate(hy_cvol)
           end if
           call deallocate_fxscr()
           
           call tileDesc%releaseDataPtr(Uin,CENTER)
           call itor%next()
        end do
     !        call Timers_stop("flux correct")
        call Grid_releaseTileIterator(itor)
     end if !Flux correction
  
  else !flux correct per level
     print *, "Flux correct per level"
     !----------------------------------------!
     !*****Flux correction per level Occurs***!
     !----------------------------------------!
     do lev=maxLev,1,-1
        
        !Once the finest level is completed, place averaged fine fluxes into
        !current coarse semipermanent flux storage (SPFS)
        if (lev < maxLev) call Grid_communicateFluxes(ALLDIR,lev)
        ! Loop over blocks and compute Hydro update block-by-block
        !~ For now tiling is disabled until we can confirm block registers are the same as tile registers
        nullify(Uin)
        call Grid_getTileIterator(itor,LEAF,level=lev,tiling=.FALSE.)
        do while(itor%isValid())
           call itor%currentTile(tileDesc)
           blkLimits(:,:)=tileDesc%limits
           blkLimitsGC(:,:)=tileDesc%blkLimitsGC
           grownLimits(:,:)=tileDesc%grownLimits
           call tileDesc%deltas(deltas)
           level=tileDesc%level
           call tileDesc%getDataPtr(Uin, CENTER)
           xLo = blkLimits(LOW,IAXIS); xHi = blkLimits(HIGH,IAXIS)
           yLo = blkLimits(LOW,JAXIS); yHi = blkLimits(HIGH,JAXIS)
           zLo = blkLimits(LOW,KAXIS); zHi = blkLimits(HIGH,KAXIS)
           
           
           hy_del=deltas
           call setLims(NGUARD-1,blkLimits,limits)
           ! DivB will technically be lagged by 1 step, but we need ghost zones to
           ! compute the gradients. I ain't doing more communication for a diagnostic...
           call allocate_scr(blkLimits,blkLimitsGC)
           call calcDivB(Uin,hy_del,blkLimits)
           call shockDetect(Uin,limits,blkLimitsGC)
           ! Allocate storage of fluxes, flux buffers, & gravity info
           call Timers_start("scratch")
           ! U* = U0
           call saveState(Uin,blkLimits,blkLimitsGC)
           call Timers_stop("scratch")
           
           !Loop stages
           do stage=1,last_stage
              ! calculate gravitational acceleration based on current value of GPOT_VAR
              ! This is stored in module-scope variable hy_grav
              
              call  hy_rk_getGraveAccel (hy_del,limits,blkLimitsGC)
              
              !Set needed number of guard cells to update based on
              !current stage for telescoping update
              call setLims(limits_array(stage), blkLimits, limits)
              
              ! Perform reconstruction and flux calculation
              ! In Stage 1, compute low-side fluxes and update for NSTENCIL guardcells
              call Timers_start("getFaceFlux")
              
              
              call hy_rk_getFaceFlux (blklimits,blkLimitsGC, limits)
              
              
              
              call Timers_stop("getFaceFlux")
              
              !In the last stage, modify fluxes on the coarse side of fine coarse boundaries.
              !This incorporates fluxes calculated in the last stage & the 'flux difference'
              !introduced on fine coarse boundaries.
              
              if (stage == last_stage) then
                 if (lev < maxLev) then
                    hy_flx(1:NFLUXES,&
                         blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                         blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                         blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))=>hya_flx
                    
                    hy_fly(1:NFLUXES,&
                         blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                         blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                         blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))=>hya_fly
                    
                    hy_flz(1:NFLUXES,&
                         blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                         blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                         blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))=>hya_flz
                    
                    call Grid_correctFluxData_xtra(tileDesc,1/weights(stage),&
                         hy_flx(:,xLo:xHi+1,yLo:yHi    ,zLo:zHi    ),&
                         hy_fly(:,xLo:xHi  ,yLo:yHi+K2D,zLo:zHi    ),&
                         hy_flz(:,xLo:xHi  ,yLo:yHi    ,zLo:zHi+K3D),&
                         blkLimits(LOW,:),-1/weights(stage),&
                         hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ)
                    nullify(hy_flx);nullify(hy_fly);nullify(hy_flz)
                 endif
              endif
              
              call addFluxes(lev,blkLimits,weights(stage),addFlux_array(stage)) 
              
              ! Now update solution based on conservative fluxes
              ! See select_RK_scheme() for explicit outline of C1, C2, & C3
              ! U* = C1 * U0 + C2 * U* + C3 * dt*L(U*)
              coeffs = coeff_array(stage,:)
              call Timers_start("updateSoln")
              if (hy_geometry /= CARTESIAN) then
                 allocate(hy_farea(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC))
                 call Grid_getCellFaceAreas(IAXIS,level,loGC,hiGC,hy_farea)
                 allocate(hy_cvol(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC))
                 call Grid_getCellVolumes(level,loGC,hiGC,hy_cvol)
                 allocate(hy_xCenter(xLoGC:xHiGC))
                 call Grid_getCellCoords(IAXIS, CENTER, level, loGC, hiGC, hy_xCenter)
                 allocate(hy_xLeft(xLoGC:xHiGC))
                 call Grid_getCellCoords(IAXIS, LEFT_EDGE, level, loGC, hiGC, hy_xLeft)
                 allocate(hy_xRight(xLoGC:xHiGC))
                 call Grid_getCellCoords(IAXIS, RIGHT_EDGE, level, loGC, hiGC, hy_xRight)
                 allocate(hy_yCenter(yLoGC:yHiGC))
                 call Grid_getCellCoords(JAXIS, CENTER, level, loGC, hiGC, hy_yCenter)
                 allocate(hy_zCenter(zLoGC:zHiGC))
                 call Grid_getCellCoords(KAXIS, CENTER, level, loGC, hiGC, hy_zCenter)
              endif
              
              call hy_rk_updateSoln(Uin,blkLimits,blklimitsGC,level,hy_del, dt, dtOld, limits, coeffs)
!!$              call tileDesc%releaseDataPtr(Uin,CENTER)
              if (hy_geometry /= CARTESIAN) then
                 deallocate(hy_xCenter)
                 deallocate(hy_xLeft)
                 deallocate(hy_xRight)
                 deallocate(hy_farea)
                 deallocate(hy_cvol)
                 deallocate(hy_yCenter)
                 deallocate(hy_zCenter)
              end if
              
              ! -------------------------------------- Deal with below here later ---------------------------------------------------__!
              
#if NSPECIES>0
              call Driver_abort("Grid_renormAbundance not implemented in SPARK with GPU offloading yet")
#endif 
              
              call Timers_stop("updateSoln")
              ! Update EOS based on intermediate solution
              
              call Timers_start("eos")
              
              call Eos_wrapped(MODE_DENS_EI,limits,hy_starState)
              
              call Timers_stop("eos")
              
              if (stage == last_stage) then
                 ! Finally, store the output and free up the scratch array
                 call Timers_start("scratch")
                 call updateState(Uin,blkLimits,blkLimitsGC)
                 call Timers_stop("scratch")
              endif
           end do!stage loop
           !Put flux buffer information into SPFS
           if (lev > 1) call Grid_putFluxData_block(tileDesc,hy_fluxBufX,hy_fluxBufY,hy_fluxBufZ,&
                blkLimits(LOW,:))
           call deallocate_scr()
           call tileDesc%releaseDataPtr(Uin,CENTER)
           call itor%next()
        end do !!block loop
        call Grid_releaseTileIterator(itor)        
     enddo!loop over levels
  endif !Flux correct per Level

! Reset local maximum hyperbolic speed. This will be updated in Hydro_computeDt.
  hy_lChyp = TINY(1.0)
  
  call Timers_stop("Hydro")
  
contains
  
#include "Hydro_funcs.F90"
  
  
  subroutine check_if_on_GPU()
!$  use omp_lib, ONLY : omp_is_initial_device
    use Driver_interface, ONLY : Driver_abort
    implicit none
    logical :: onCPU
    
    onCPU = .TRUE.
    
    !$omp target map(tofrom:onCPU)
    !$  onCPU = omp_is_initial_device()
    !$omp end target
    
    if (onCPU) then
       print *, "---------------- Running on CPU --------------------------------"
       ! call Driver_abort("Unable to run on GPU")
    else
       print *, "---------------- Running on GPU --------------------------------"
    end if
    
  end subroutine check_if_on_GPU
  
end subroutine Hydro
