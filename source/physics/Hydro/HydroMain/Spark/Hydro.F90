!!****if* source/physics/Hydro/HydroMain/Spark/Hydro
!!
!!
!! NAME
!!
!!  Hydro
!!
!!
!! SYNOPSIS
!!
!!  Hydro(real(IN)    :: timeEndAdv,
!!        real(IN)    :: dt,
!!        real(IN)    :: dtOld,
!!        integer(IN) :: sweepOrder)
!!
!!
!! DESCRIPTION
!!
!!  Performs physics update in a directionally unsplit fashion.
!!
!!  dt gives the timestep through which this update should advance,
!!  and timeEndAdv tells the time that this update will reach when
!!  it finishes.  dtOld gives the previously taken timestep.
!!
!! ARGUMENTS
!!
!!  timeEndAdv - end time
!!  dt         - timestep
!!  dtOld      - old timestep
!!  sweepOrder - dummy argument for the unsplit scheme, just a dummy
!!               variable to be consistent with a toplayer stub function
!!
!!***
!!Reorder(4): hy_fl[xyz]

subroutine Hydro(timeEndAdv, dt, dtOld, sweepOrder)
  use Hydro_data, ONLY : hy_useHydro, hy_fluxCorrect, hy_fluxCorrectPerLevel, hy_starState, &
      hy_gcMask, hy_lChyp, hy_C_hyp, hy_globalComm, &
      hy_flx, hy_fly, hy_flz, hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ,hy_tiny, hy_hybridRiemann
  use hy_rk_interface, ONLY : hy_rk_eos, hy_rk_getFaceFlux, hy_rk_getGravAccel, hy_rk_updateSoln
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_communicateFluxes, &
      Grid_fillGuardCells, Grid_getTileIterator, &
      Grid_releaseTileIterator, Grid_getMaxRefinement, &
      Grid_correctFluxData_xtra, Grid_putFluxData_block, &
      Grid_getFluxCorrData_block
  use Eos_interface, ONLY : Eos_wrapped
  use IO_interface, ONLY : IO_setScalar
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t
 
#include "Simulation.h"
#include "constants.h"
#include "Flashx_mpi_implicitNone.fh"

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: blockDesc
  real,    intent(in) :: timeEndAdv, dt, dtOld
  integer, intent(IN) :: sweepOrder
  integer :: n, error, level, maxLev=-1
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: xLo,xHi,yLo,yHi,zLo,zHi
  integer, dimension(LOW:HIGH,MDIM) :: limits
  real :: hdt
  integer :: stage, last_stage
  real, dimension(3) :: coeffs, weights
  integer, dimension(3) :: limits_array
  real, dimension(3,3) :: coeff_array
  logical, dimension (3) :: addFlux_array
  integer :: i,j,k,v

  if (.NOT. hy_useHydro) return

  call Timers_start("Hydro")

#ifdef OMP_OL
  call check_if_on_GPU()
#endif

  hdt = 0.5*dt

  !Array indicating whether to add flux into flux buffers (True)
  !or overwrite them (False).
  !This array will work for RK2 & RK3
  addFlux_array = (/.false.,.true.,.true./)
 
  !set up quantities specific to RK scheme (lives in Hydro_funcs)
  call select_RK_scheme(coeff_array,last_stage,limits_array,weights)

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

!!!*** here we include keyword about perlevel or all levels at once

!!!***from here code snippet starts
  call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.false.,maskSize=NUNK_VARS,mask=hy_gcMask)
  !-------------------------------------------------------------------!
  !***NO Flux correction    or   Flux correction but NOT per level****!
  !-------------------------------------------------------------------!
  if ((.NOT.hy_fluxCorrect).OR.((hy_fluxCorrect).AND.(.NOT.hy_fluxCorrectPerLevel))) then
!!!*** keyword to generate information for the runtime orchestrator, indicate bookend
    ! Loop over blocks and compute Hydro update block-by-block

    call Grid_getTileIterator(itor,LEAF,tiling=.FALSE.) !FALSE? Try True?
    do while(itor%isValid())
       call itor%currentTile(blockDesc)

       ! DivB will technically be lagged by 1 step, but we need ghost zones to
       ! compute the gradients. I ain't doing more communication for a diagnostic...
       call setLims(limits, NGUARD-1,blockDesc)
!!!*** this is where the first part of code snippet ends
!!!*** for all of these functions we decide if we should replace passing blockDesc by arguments     

       !---JaredC--- Action 2 (no dependence within interation) 
       call calcDivB(blockDesc) !---JaredC--- This could be done on GPU or CPU and asynchronous to the rest of the sections
                                !---JaredC--- This could be pulled out of the loop and given its own loop.


!!!*** there is some scratch data need, so may need composer help
      !  call shockDetect(blockDesc,limits)
       ! Setup scratch storage of block data (grav included)
       !---JaredC--- Send data to device and/or allocate arrays
       call Timers_start("Offloaded Section")
       call Timers_start("scratch")
       ! U* = U0
      !  The data is sent to the GPU at this point.
       call saveState(blockDesc)
       call Timers_stop("scratch")

       !---JaredC--- Action 3 
       !Begin loop over stages
       do stage=1,last_stage
         ! calculate gravitational acceleration based on current value of GPOT_VAR
         ! This is stored in module-scope variable hy_grav
!!!*** only needs transpilation
         call hy_rk_getGravAccel(blockDesc,limits)

         !Set needed number of guard cells to update based on
         !current stage for telescoping update
         call setLims(limits,limits_array(stage), blockDesc)
       
         call Timers_start("getFaceFlux")
!!!*** this functions has some allocation and data rearrangement so need intervention from composer
          call hy_rk_getFaceFlux(blockDesc, limits)

         call Timers_stop("getFaceFlux")
!------- Add this in -----
        !  if (hy_fluxCorrect) call addFluxes(weights(stage),addFlux_array(stage)) 
         
         ! Now update solution based on conservative fluxes
         ! See select_RK_scheme() for explicit outline of C1, C2, & C3
         ! U* = C1 * U0 + C2 * U* + C3 * dt*L(U*)

         !Set proper coefficients for the given stage
         coeffs = coeff_array(stage,:)

         call Timers_start("updateSoln")
!!!*** just needs transpiler
         call hy_rk_updateSoln(blockDesc,dt,dtOld,limits,coeffs)
         call Timers_stop("updateSoln")
         ! Update EOS based on intermediate solution
         call Timers_start("eos")
!!!*** There is invocation to Eos_get/putData so needs composer help     
#ifdef OMP_OL
         call hy_rk_eos_offloaded(limits)
#else
         call hy_rk_eos(limits)
#endif
         call Timers_stop("eos")

         !Finally update the state
         if (stage == last_stage) then
          call Timers_start("scratch")
          ! The data is sent off the GPU at this point
          call updateState(blockDesc)
          call Timers_stop("scratch")
         endif

       enddo!stage loop
       !Store flux buffer in semipermanent flux storage (SPFS) 
      !  if (hy_fluxCorrect) call Grid_putFluxData_block(blockDesc,&
      !                           hy_fluxBufX,hy_fluxBufY,hy_fluxBufZ,blkLimits(LOW,:)) 
       call Timers_stop("Offloaded Section")
       call itor%next()
    enddo!block loop
    call Grid_releaseTileIterator(itor)   
!*** keyword to generate information for the runtime orchestrator, indicate bookend

!     if (hy_fluxCorrect) then
!        !Communicate the fine fluxes
!        call Grid_communicateFluxes(ALLDIR,UNSPEC_LEVEL)
 
!        call Timers_start("flux correct")
!        ! Loop over blocks and correct block-edge solution
! !!!*** keyword to generate information for the runtime orchestrator, indicate bookend
!        call Grid_getTileIterator(itor,LEAF,tiling=.FALSE.) !FALSE? Try True
!        do while(itor%isValid())
!           call itor%currentTile(blockDesc)
       
!           blkLimits(:,:)   = blockDesc%limits
!           blkLimitsGC(:,:) = blockDesc%blkLimitsGC
!           blkLimits(:,1:NDIM) = blkLimits(:,1:NDIM)
!           blkLimitsGC(:,1:NDIM) = blkLimitsGC(:,1:NDIM) 
!           !Get 'Flux deltas' on coarse side of fine coarse boundaries; 
!           !all other values are 0.
!           call Grid_getFluxCorrData_block(blockDesc,hy_fluxBufX,hy_fluxBufY,hy_fluxBufZ,&
!                                           blkLimits(LOW,:))
! !!!*** only needs transpiler
!           call hy_rk_correctFluxes(blockDesc,dt)
!           call itor%next()
!        end do
!        call Timers_stop("flux correct")
!        call Grid_releaseTileIterator(itor)
!     end if !Flux correction

  else !flux correct per level
    print *, "Flux correct per level"
    !----------------------------------------!
    !*****Flux correction per level Occurs***!
    !----------------------------------------!
    do level=maxLev,1,-1
!!!*** keyword to generate information for the runtime orchestrator, indicate bookend
      !Once the finest level is completed, place averaged fine fluxes into
      !current coarse semipermanent flux storage (SPFS)
      if (level < maxLev) call Grid_communicateFluxes(ALLDIR,level)
      ! Loop over blocks and compute Hydro update block-by-block
      !~ For now tiling is disabled until we can confirm block registers are the same as tile registers
      call Grid_getTileIterator(itor,LEAF,level=level,tiling=.FALSE.)
      do while(itor%isValid())
         call itor%currentTile(blockDesc)
         blkLimits(:,:)   = blockDesc%limits
         blkLimitsGC(:,:) = blockDesc%blkLimitsGC 
         blkLimits(:,1:NDIM) = blkLimits(:,1:NDIM)
         blkLimitsGC(:,1:NDIM) = blkLimitsGC(:,1:NDIM)
         xLo = blkLimits(LOW,IAXIS); xHi = blkLimits(HIGH,IAXIS)
         yLo = blkLimits(LOW,JAXIS); yHi = blkLimits(HIGH,JAXIS)
         zLo = blkLimits(LOW,KAXIS); zHi = blkLimits(HIGH,KAXIS)

         call setLims(limits, NGUARD-1, blockDesc)
!!!*** this is where the first part of code snippet ends
!!!*** for all of these functions we decide if we should replace passing blockDesc by arguments     
         ! DivB will technically be lagged by 1 step, but we need ghost zones to
         ! compute the gradients. I ain't doing more communication for a diagnostic...         
         call calcDivB(blockDesc)

!!!*** there is some scratch data need, so may need composer help
         call shockDetect(blockDesc,limits)
         ! Allocate storage of fluxes, flux buffers, & gravity info
         call Timers_start("scratch")
         ! U* = U0
         call saveState(blockDesc)
         call Timers_stop("scratch")
        
         !Loop stages
         do stage=1,last_stage
           ! calculate gravitational acceleration based on current value of GPOT_VAR
           ! This is stored in module-scope variable hy_grav
!!!*** only needs transpilation
           call hy_rk_getGravAccel(blockDesc,limits)

           !Set needed number of guard cells to update based on
           !current stage for telescoping update
           call setLims(limits,limits_array(stage),blockDesc)
           
           ! Perform reconstruction and flux calculation
           ! In Stage 1, compute low-side fluxes and update for NSTENCIL guardcells
           call Timers_start("getFaceFlux")
!!!*** this functions has some allocation and data rearrangement so need intervention from composer
           call hy_rk_getFaceFlux(blockDesc, limits)
           call Timers_stop("getFaceFlux")

           !In the last stage, modify fluxes on the coarse side of fine coarse boundaries.
           !This incorporates fluxes calculated in the last stage & the 'flux difference'
           !introduced on fine coarse boundaries.

           if (stage == last_stage) then
             if (level < maxLev) then
               call Grid_correctFluxData_xtra(blockDesc,1/weights(stage),&
                    hy_flx(:,xLo:xHi+1,yLo:yHi    ,zLo:zHi    ),&
                    hy_fly(:,xLo:xHi  ,yLo:yHi+K2D,zLo:zHi    ),&
                    hy_flz(:,xLo:xHi  ,yLo:yHi    ,zLo:zHi+K3D),&
                    blkLimits(LOW,:),-1/weights(stage),&
                    hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ)
             endif
           endif

           call addFluxes(weights(stage),addFlux_array(stage)) 

           ! Now update solution based on conservative fluxes
           ! See select_RK_scheme() for explicit outline of C1, C2, & C3
           ! U* = C1 * U0 + C2 * U* + C3 * dt*L(U*)
           coeffs = coeff_array(stage,:)
           call Timers_start("updateSoln")

!!!*** just needs transpiler
           call hy_rk_updateSoln(blockDesc,dt,dtOld,limits,coeffs)
           call Timers_stop("updateSoln")
           ! Update EOS based on intermediate solution
        
           call Timers_start("eos")
!!!*** There is invocation to Eos_get/putData so needs composer help     
           call hy_rk_eos(limits)
           call Timers_stop("eos")

           if (stage == last_stage) then
             ! Finally, store the output and free up the scratch array
             call Timers_start("scratch")
             call updateState(blockDesc)
             call Timers_stop("scratch")
           endif
         end do!stage loop
         !Put flux buffer information into SPFS
         if (level > 1) call Grid_putFluxData_block(blockDesc,hy_fluxBufX,hy_fluxBufY,hy_fluxBufZ,&
                              blkLimits(LOW,:))
         
         call itor%next()
      end do!block loop
      call Grid_releaseTileIterator(itor)
!!!*** keyword to generate information for the runtime orchestrator, indicate bookend
    enddo!loop over levels
  endif !Flux correct per Level
! !!!###
!!!*** keyword to generate information for the runtime orchestrator, indicate bookend  
!!!*** this is where the second part of the composer handing per level or not comes in

  ! Reset local maximum hyperbolic speed. This will be updated in Hydro_computeDt.
  hy_lChyp = TINY(1.0)

  call Timers_stop("Hydro")

contains

#include "Hydro_funcs.F90"


subroutine check_if_on_GPU()
  use omp_lib, ONLY : omp_is_initial_device
  use Driver_interface, ONLY : Driver_abort
  implicit none
  logical :: onCPU

  onCPU = .TRUE.

  !$omp target map(tofrom:onCPU)
 !! onCPU = omp_is_initial_device()
  !$omp end target
  
  if (onCPU) then
      print *, "---------------- Running on CPU --------------------------------"
      ! call Driver_abort("Unable to run on GPU")
  else
    print *, "---------------- Running on GPU --------------------------------"
  end if
  
  end subroutine check_if_on_GPU

end subroutine Hydro
