!  Apply parameterized, multi-stage progress variable burner 
!  to all blocks in specified list.
!
!   Dean Townsley 2008
!

subroutine Burn ( blockCount, blockList, dt )
  
  use Grid_interface, ONLY  : Grid_fillGuardCells, &
       Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getSingleCellVol
  use Eos_interface, ONLY   : Eos_wrapped
  use Hydro_interface, ONLY : Hydro_detectShock
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use bn_paraInterface, ONLY : bn_paraBurn, bn_paraSpark, bn_paraAllSpark, &
                               bn_paraAllIgnite
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Driver_interface, ONLY : Driver_abort, Driver_getSimTime
  use Flame_interface, ONLY : Flame_getWidth

  use Burn_data, ONLY : bn_meshMe, bn_useBurn, bn_useShockBurn, bn_thermalReact, &
                        bn_autoDDT, pbIgnRad, pbIgnPhfa, pbIgnRho, pbIgnNum, &
                        pbIgnTime, pbIgnX, pbIgnY, pbIgnZ, pbIgnR, &
                        pbIgnNumMax, bn_neutLossThisProcStep, &
                        burn_lun, bn_detIgnFileName, &
                        bn_thermalReactInFlameThreshold
  
  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"

  integer, INTENT(in)                        :: blockCount
  integer, INTENT(in), DIMENSION(blockCount) :: blockList
  real,    INTENT(in)                        :: dt

  integer                    :: i, j, k, iref, l, n
  integer                    :: blockID, istat
  real                       :: eint, flame, flamedot, qdot, edotnu, time, dx

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, allocatable, dimension(:)         :: xCoord, yCoord, zCoord
  integer                                 :: xSizeCoord, ySizeCoord, zSizeCoord
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: shock
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: react_proximity
  real, pointer, dimension(:,:,:,:)                    :: solnData

  real :: flamewidth
  integer :: ii, jj, kk, maxdi

  logical :: ignition_test, bcastMask, spark_block
  real :: spark_block_radius, radius
  logical, allocatable, dimension(:) :: ignition_conditions
  integer :: ignition_num
  real, allocatable, dimension(:) :: ignition_coords

  integer, allocatable, dimension(:) :: broadcast_order
  integer :: det_num
  real, allocatable, dimension(:) :: det_xCoord, det_yCoord, det_zCoord
  real :: detX, detY, detZ, dist, phi_fa_det

  logical :: ignite_detonation

  logical :: gcMask(NUNK_VARS)

#ifndef EINT_VAR
  real :: energyKinetic
#endif

  integer :: point(MDIM)
  real    :: dvol

100 format("*** igniting at ",f7.5,4(1x,es10.3)," ***")
  ! --------------------------------------------------------
  ! 0. do nothing if burn is turned off and initalize some things
  if (.not. bn_useBurn) return
  
  call Timers_start("burn")
  
  call Driver_getSimTime(time)
  call Flame_getWidth(flamewidth)


  ! ---------------------------------------------------
  ! 1. bring guardcells up to date
  ! ---------------------------------------------------
  gcMask = .FALSE.
  ! If burning in shocks is disabled (the normal case)
  ! we will need the pressure and velocity in order to do shock detection
  if (bn_thermalReact .and. (.NOT. bn_useShockBurn)) then
     gcMask(PRES_VAR) = .TRUE.
     gcMask(VELX_VAR) = .TRUE.
#if NDIM > 1
     gcMask(VELY_VAR) = .TRUE.
#endif
#if NDIM > 2
     gcMask(VELZ_VAR) = .TRUE.
#endif
  endif
  
  if (bn_thermalReact) then
     ! need guardcells for burning proximity checks
     gcMask(PHFA_MSCALAR) = .true.

     ! fill guardcells
     ! we don't need to do the EOS because if there is a shock at the
     ! edge of the block, then the next block should be fully refined also and
     ! thus at the same refinement as this one.  The EOS call is only
     ! really needed when the neighboring block is at a different
     ! refinement level.
     call Grid_fillGuardCells(CENTER, ALLDIR, &
       doEos=.false., maskSize=NUNK_VARS,  mask = gcMask, makemaskConsistent=.true.)
  endif

  if (bn_autoDDT .and. bn_thermalReact) then

     ! -------------------------------------------------------
     ! 2.  iterate over blocks and check for ignition conditions 
     ! -------------------------------------------------------
     ! determine the maximum number of ignition points to allow

     ! initialize
     ignition_num = 0
     allocate(ignition_coords(NDIM*pbIgnNumMax),stat=istat)
     if (istat/=0) call Driver_abortFlash("Unable to allocate ignition_coords")
  
     ! store the number of ignition points in index zero
     do n = 1, blockCount

        spark_block = .false.
        spark_block_radius = 0.e0
        blockID = blockList(n)  

        ! -------------------------------
        ! 2.1 initialize quantities for this block
        ! -------------------------------
        call Grid_getBlkPtr(blockID,solnData)
     
        ! get coordinate positions
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        xSizeCoord = blkLimitsGC(HIGH,IAXIS)
        ySizeCoord = blkLimitsGC(HIGH,JAXIS)
        zSizeCoord = blkLimitsGC(HIGH,KAXIS)
        allocate(xCoord(xSizeCoord))
        allocate(yCoord(ySizeCoord))
        allocate(zCoord(zSizeCoord))
        call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,xCoord,xSizeCoord)
        call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,yCoord,ySizeCoord)
        call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,zCoord,zSizeCoord)

        ! --------------------------------
        ! 2.2 loop over all interior zones and check for ignition conditions 
        ! --------------------------------
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
                 ! --------------------------------
                 ! 2.2.1 check for ignition conditions
                 ! --------------------------------
#ifndef FLAM_MSCALAR
                 flame   = 0.e0
#else
                 flame   = solnData(FLAM_MSCALAR,i,j,k)
#endif
                 call bn_paraSpark( xCoord(i), yCoord(j), zCoord(k), &
                                    solnData(DENS_VAR,i,j,k),        &
                                    solnData(PRES_VAR,i,j,k),        &
                                    solnData(PHFA_MSCALAR,i,j,k),    &
                                    flame,                           &
                                    solnData(CI_MSCALAR,i,j,k),      &
                                    solnData(NEI_MSCALAR,i,j,k),     &
                                    detX, detY, detZ, ignition_test )
                 ! --------------------------------
                 ! 2.2.2 store ignition coordinates
                 ! --------------------------------
                 if ( ignition_test ) then

                    radius = sqrt( detX**2 + detY**2 + detZ**2 )

                    ! check to see if we have already ignited on this block
                    if (spark_block) then

                       ! check which is further out from the center
                       if (radius > spark_block_radius) then

                          ! if current radius is larger, then overwrite previous
                          ! ignition point
                          spark_block_radius = radius
                          ignition_coords((ignition_num-1)*NDIM + 1) = detX
#if NDIM >= 2
                          ignition_coords((ignition_num-1)*NDIM + 2) = detY
#endif
#if NDIM > 2
                          ignition_coords((ignition_num-1)*NDIM + 3) = detZ
#endif
                       endif !! (radius > spark_block_radius)

                       ! otherwise, our previous ignition point is further out
                       ! and we will keep it

                    else !! (spark_block)

                       ! we have not already found an ignition condition on this
                       ! block so, lets set a new one

                       ! make sure we have enough space
                       if ( ignition_num < pbIgnNumMax ) then

                          spark_block = .true.
                          ignition_coords( ignition_num*NDIM + 1 ) = detX
#if NDIM >= 2
                          ignition_coords( ignition_num*NDIM + 2 ) = detY
#endif
#if NDIM > 2
                          ignition_coords( ignition_num*NDIM + 3 ) = detZ
#endif
                          ignition_num = ignition_num + 1

                       ! we have too many ignition points, 
                       ! need to create more space
                       else !! ( ignition_num < pbIgnNumMax )

                         call Driver_abortFlash("Not enough space to store all ignition points")

                       endif !! ( ignition_num < pbIgnNumMax )

                    endif !! (spark_block)

                 endif !! (ignition_test)

              enddo
           enddo
        enddo

        ! --------------------------------
        ! 2.3 Finish up:
        !     Release stuff
        ! --------------------------------
        call Grid_releaseBlkPtr(blockID,solnData)
        deallocate(xCoord)
        deallocate(yCoord)
        deallocate(zCoord)
     enddo

     ! -------------------------------------------------------
     ! 3.  communicate ignition points
     ! -------------------------------------------------------
     ! store number of ignition points in first index
     if (ignition_num .eq. 0) then
        deallocate(ignition_coords)
     endif

     ! allocate space for det coords
!     allocate(det_xCoord(pbIgnNumMax),stat=istat)
!     if (istat/=0) call Driver_abortFlash("Unable to allocate det_xCoord")
!     allocate(det_yCoord(pbIgnNumMax),stat=istat)
!     if (istat/=0) call Driver_abortFlash("Unable to allocate det_yCoord")
!     allocate(det_zCoord(pbIgnNumMax),stat=istat)
!     if (istat/=0) call Driver_abortFlash("Unable to allocate det_zCoord")

     ! communicate detonation points and consolidate
     ! this means only det points which are not near each other and
     ! not near previous det points are returned, we still have to 
     ! check these det points for being in ash
     det_num = ignition_num
     call bn_paraAllSpark( ignition_coords, det_num, &
                           det_xCoord, det_yCoord, det_zCoord )

     ! bn_paraAllSpark has deallocated ignition_coords and has
     ! allocated det_xCoord, det_yCoord, det_zCoord

     ! -------------------------------------------------------
     ! 4.  iterate over detonation points and blocks to
     !     see if detonation points are valid
     ! -------------------------------------------------------
     ! make sure there are detonation points
     if (det_num > 0) then

        allocate(ignition_conditions(det_num),stat=istat)
        if (istat/=0)  &
           call Driver_abortFlash("Unable to allocate ignition_conditions")
        ignition_conditions(:) = .true.
 
        det_search:  do l = 1, det_num
           ! ---------------------------
           ! 4.1 iterate over all blocks
           ! ---------------------------
           do n = 1, blockCount
    
              blockID = blockList(n)

              ! -------------------------------
              ! 4.1.1 initialize quantities for this block
              ! -------------------------------
              call Grid_getBlkPtr(blockID,solnData)
     
              ! get coordinate positions, used for shock detection and
              ! non-flame burning proximity detection
              call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
              xSizeCoord = blkLimitsGC(HIGH,IAXIS)
              ySizeCoord = blkLimitsGC(HIGH,JAXIS)
              zSizeCoord = blkLimitsGC(HIGH,KAXIS)
              allocate(xCoord(xSizeCoord))
              allocate(yCoord(ySizeCoord))
              allocate(zCoord(zSizeCoord))
              call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,xCoord,xSizeCoord)
              call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,yCoord,ySizeCoord)
              call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,zCoord,zSizeCoord)

              ! --------------------------------
              ! 4.1.2 loop over all interior zones
              ! --------------------------------
              do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                 do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                       ! --------------------------------
                       ! 4.1.2.1 check if we are near a detonation point
                       ! --------------------------------
                       dist = ( det_xCoord(l) - xCoord(i) )**2
#if NDIM >= 2
                       dist = dist + ( det_yCoord(l) - yCoord(j) )**2
#endif
#if NDIM > 2
                       dist = dist + ( det_zCoord(l) - zCoord(k) )**2
#endif
                       dist = sqrt( dist )

                       if (dist <= pbIgnRad) then

                          ! check if we are in ash, and we are at ignition density
                          if ( solnData(PHFA_MSCALAR,i,j,k) > pbIgnPhfa ) then
   
                             ! we are in ash 
                             ! so invalidate detonation point
                             ignition_conditions(l) = .false.

                             ! we are done so finish up with this block
                             call Grid_releaseBlkPtr(blockID,solnData)
                             deallocate(xCoord)
                             deallocate(yCoord)
                             deallocate(zCoord)

                             ! move to the next detonation point
                             cycle det_search
                          endif
                       endif

                    enddo
                 enddo
              enddo
     
              ! --------------------------------
              ! 4.1.3 Finish up:
              !     Release stuff
              ! --------------------------------
              call Grid_releaseBlkPtr(blockID,solnData)
              deallocate(xCoord)
              deallocate(yCoord)
              deallocate(zCoord)
           enddo

        enddo det_search

        ! -----------------------
        ! 4.2 communicate any invalidated detonation points
        ! -----------------------
        call bn_paraAllIgnite( ignition_conditions, det_num )

        ! announce detonation points and save to global list
        if (pbIgnNum + count(ignition_conditions) < pbIgnNumMax) then
           do l = 1, det_num
              if ( ignition_conditions(l) ) then
                 if (bn_meshMe .eq. MASTER_PE) then
                    detX = det_xCoord(l)
#if NDIM >= 2
                    detY = det_yCoord(l)
#else
                    detY = 0.e0
#endif
#if NDIM > 2
                    detZ = det_zCoord(l)
#else
                    detZ = 0.e0
#endif
                    open (burn_lun, file=bn_detIgnFileName, &
                          position='append', iostat=istat)
                    if (istat/=0) call Driver_abortFlash("Unable to open detonation ignition points file")

                    write (burn_lun,*) time, detX, detY, detZ

                    close(unit=burn_lun)

                 endif

                 pbIgnNum = pbIgnNum + 1
                 pbIgnTime( pbIgnNum ) = time
                 pbIgnX( pbIgnNum ) = det_xCoord(l)
#if NDIM >= 2
                 pbIgnY( pbIgnNum ) = det_yCoord(l)
#else
                 pbIgnY( pbIgnNum ) = 0.e0
#endif
#if NDIM > 2
                 pbIgnZ( pbIgnNum ) = det_zCoord(l)
#else
                 pbIgnZ( pbIgnNum ) = 0.e0
#endif
              endif
           enddo
        else
           call Driver_abortFlash("Not enough space to save detonation points")
        endif

     endif

  endif

  ! ------------------------------------------------------
  ! 5.  iterate over blocks, apply burning evolution and update quantities
  ! -------------------------------------------------------
  do n = 1, blockCount
     
     blockID = blockList(n)

     ! -------------------------------
     ! 5.1 initialize quantities for this block
     !     including shock and proximity, which require unchanged neighbor cells
     ! -------------------------------
     call Grid_getBlkPtr(blockID,solnData)
     
     ! get coordinate positions, used for shock detection and
     ! non-flame burning proximity detection
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     xSizeCoord = blkLimitsGC(HIGH,IAXIS)
     ySizeCoord = blkLimitsGC(HIGH,JAXIS)
     zSizeCoord = blkLimitsGC(HIGH,KAXIS)
     allocate(xCoord(xSizeCoord))
     allocate(yCoord(ySizeCoord))
     allocate(zCoord(zSizeCoord))
     call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,xCoord,xSizeCoord)
     call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,yCoord,ySizeCoord)
     call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,zCoord,zSizeCoord)

     dx = xCoord(2)-xCoord(1) ! also assume square grid
     
     ! shock detect if burning is turned off in shocks
     if (bn_thermalReact .and. (.NOT. bn_useShockBurn)) then
        call Hydro_detectShock(solnData, shock, blkLimits, blkLimitsGC, (/0,0,0/), &
             xCoord,yCoord,zCoord)
     else
        shock(:,:,:) = 0
     endif

#ifdef SHK_VAR
     solnData(SHK_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
                      blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                      blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
                shock(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
                      blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),   &
                      blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
#endif
     
     !  Check for proximity of a reacting region for each cell
     !  this is used to help control thermal burning inside flame
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              !  proximity is in units of flame width (rounded up to nearest cell)
              ! no need to sqrt react_proximity because we are comparing to 1.0
              react_proximity(i,j,k) = 2.0  ! > 1 supresses reaction in flame
! allow burning with no flame module compiled in
! proximity doesn't matter with no flame, but be safe and just set leave it default
#ifdef FLAM_MSCALAR
              if (bn_thermalReact) then
                 maxdi = int(ceiling(flamewidth/dx))
                 do ii = -maxdi, maxdi
                    do jj = -maxdi*K2D, maxdi*K2D
                       do kk = -maxdi*K3D, maxdi*K3D
                          if ( solnData(PHFA_MSCALAR,i+ii,j+jj,k+kk) - solnData(FLAM_MSCALAR,i+ii,j+jj,k+kk) &
                                     > bn_thermalReactInFlameThreshold ) then
                             react_proximity(i,j,k) = min(react_proximity(i,j,k),real(ii**2+jj**2+kk**2)/maxdi**2)
                          endif
                       enddo
                    enddo
                 enddo
              endif
#endif
           enddo
        enddo
     enddo

     ! --------------------------------
     ! 5.2 loop over all interior zones and apply evolution
     ! --------------------------------
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              ! --------------------------------
              ! 5.2.1 set up internal energy, flame inputs
              !       and check detonation ignition points
              ! --------------------------------
! eint can be eliminated, leaving just ENER
#ifdef EINT_VAR
              eint    = solnData(EINT_VAR,i,j,k)
#else
              energyKinetic = solnData(VELX_VAR,i,j,k)**2
#if NDIM >= 2
              energyKinetic = energyKinetic + solnData(VELY_VAR,i,j,k)**2
#endif
#if NDIM > 2
              energyKinetic = energyKinetic + solnData(VELZ_VAR,i,j,k)**2
#endif
              eint    = solnData(ENER_VAR,i,j,k) - 0.5*energyKinetic
#endif
              
! allow burning with no flame module compiled in
#ifndef FLAM_MSCALAR
              flame   = 0.e0
              flamedot = 0.0
#else
              flame   = solnData(FLAM_MSCALAR,i,j,k)
              flamedot = solnData(FLDT_VAR,i,j,k)
#endif
               ! --------------------------------
              ! 5.2.2 check if we are near a detonation point
              ! --------------------------------
              ignite_detonation = .false.
              if (bn_autoDDT) radius = pbIgnRad
              do l = 1, pbIgnNum

                 dist = (pbIgnX(l) - xCoord(i))**2
#if NDIM >= 2
                 dist = dist + (pbIgnY(l) - yCoord(j))**2
#endif
#if NDIM > 2
                 dist = dist + (pbIgnZ(l) - zCoord(k))**2
#endif
                 dist = sqrt(dist)

                 if (.not. bn_autoDDT) radius = pbIgnR(l)
                 if ( pbIgnTime(l) < time+dt .and. pbIgnTime(l) >= time &
                      .and. dist <= radius ) then
                    write (6,100) pbIgnTime(l), pbIgnX(l), pbIgnY(l), pbIgnZ(l), radius
                    ignite_detonation = .true.
                    phi_fa_det = 1.0e0
                 endif
              enddo

              ! --------------------------------
              ! 2.2.2 evolve progress variables and update NSE grid quantities
              ! --------------------------------
              call bn_paraBurn( solnData(DENS_VAR,i,j,k), &
                                solnData(TEMP_VAR,i,j,k), &
                                eint, &
                                solnData(PRES_VAR,i,j,k), &
                                solnData(CI_MSCALAR,i,j,k), &
                                solnData(NEI_MSCALAR,i,j,k), &
                                flame, &
                                flamedot, &
                                solnData(PHFA_MSCALAR,i,j,k), &
                                solnData(PHAQ_MSCALAR,i,j,k), &
                                solnData(PHQN_MSCALAR,i,j,k), &
                                solnData(YE_MSCALAR,i,j,k), &
                                solnData(DYQN_MSCALAR,i,j,k), &
                                solnData(DQQN_MSCALAR,i,j,k), &
                                qdot, edotnu, dt, &
                                react_proximity(i,j,k), shock(i,j,k), &
                                ignite_detonation, &
                                phi_fa_det )
              ! save energy deposition rate
              solnData(ENUC_VAR,i,j,k)     = qdot 
              ! deposit energy
              solnData(ENER_VAR,i,j,k)     = solnData(ENER_VAR,i,j,k) + qdot*dt
#ifdef EINT_VAR
              solnData(EINT_VAR,i,j,k)     = solnData(EINT_VAR,i,j,k) + qdot*dt
#endif

              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k
              call Grid_getSingleCellVol(blockID, EXTERIOR, point, dvol) 
              bn_neutLossThisProcStep = bn_neutLossThisProcStep + solnData(DENS_VAR,i,j,k)*dvol*edotnu*dt

           enddo
        enddo
     enddo
     
     
     ! --------------------------------
     ! 2.3 Finish up:
     !     Update interior EOS quantities for this block and release stuff
     ! --------------------------------
     call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)

     call Grid_releaseBlkPtr(blockID,solnData)
     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)
  end do

  if (bn_autoDDT .and. bn_thermalReact) then

     ! deallocate stuff
     if (det_num > 0) then
        deallocate(ignition_conditions)  
        deallocate(det_xCoord)
#if NDIM >= 2
        deallocate(det_yCoord)
#endif
#if NDIM > 2
        deallocate(det_zCoord)
#endif
     endif

  endif

  call Timers_stop("burn")
  
  return
  
end subroutine Burn
