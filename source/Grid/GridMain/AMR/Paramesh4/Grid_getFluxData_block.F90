!!****if* source/Grid/GridMain/AMR/Paramesh4/Grid_getFluxData_block
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
!!  Grid_getFluxData_block
!!
!! SYNOPSIS
!!
!!  call Grid_getFluxData_block(type(Grid_tile_t)(in) :: blockDesc,
!!                              real(OUT),TARGET      :: fluxBufX(:,lo(1): ,lo(2): ,lo(3): ),
!!                              real(OUT),TARGET      :: fluxBufY(:,lo(1): ,lo(2): ,lo(3): ),
!!                              real(OUT),TARGET      :: fluxBufZ(:,lo(1): ,lo(2): ,lo(3): ),
!!                              integer(in)           :: lo(3),
!!                              integer(in)           :: axis,
!!                              logical, intent(IN), OPTIONAL  :: isFluxDensity)
!!
!! DESCRIPTION
!!
!!   Get (corrected) flux data from semipermanent flux storage (SPFS).
!!
!!    fluxBuf := "communicated fine fluxes"   AT           coarse side of f/c bdry;
!!            := "saved coarse fluxes"        ELSEWHERE AT cells touching block bdry;
!!            undef                           ELSEWHERE.
!!
!!   A stub that does nothing is used for the Uniform Grid implementation.
!!
!! ARGUMENTS
!!
!!   blockDesc : describes the current block.
!!               Note that this should be a full block, not a tile representing
!!               a partial block.
!!
!!   fluxBufX :  buffer for fluxes in IAXIS-direction
!!
!!   fluxBufY :  buffer for fluxes in JAXIS-direction; output undefined if NDIM < 2
!!
!!   fluxBufZ :  buffer for fluxes in KAXIS-direction; output undefined if NDIM < 3
!!
!!   lo :        lower bounds for the spatial indices of the flux buffers
!!
!!  axis : integer value specifying on which cell faces to get fluxes. 
!!         The options are IAXIS, JAXIS, KAXIS, or ALLDIR defined in constants.h
!!
!!   isFluxDensity : indicates, for each flux component, whether the component
!!                   is a flux proper (if TRUE) or a flux density (otherwise).
!!                   This may be either removed, or changed into a scalar flag,
!!                   later.
!!
!! NOTES
!!
!!   The arrays fluxBufX, fluxBufY, fluxBufZ are subject to index reordering.
!!
!!   flux buffer arrays should contain space for fluxes of all valid cells
!!   in the block, excluding guard cells.
!!
!!   This interface does not require level-wide fluxes to be allocated.
!!
!!   SPFS means semi-permanent flux storage. When using a Grid
!!   implementation based on AMReX, SPFS is implemented by an AMReX
!!   flux register class, such as FlashFluxRegister.
!!
!!   This implementation is specific to Paramesh 4.
!!
!! SEE ALSO
!!
!!   Grid_putFluxData_block
!!   Grid_communicateFluxes
!!   Grid_correctFluxData
!!   Hydro
!!***

!!REORDER(4): fluxBuf[XYZ]
!!REORDER(5): flux_[xyz], gr_[xyz]flx
!!REORDER(4): fluxx,fluxy,fluxz
!!REORDER(5): gr_xflx_[yz]face, gr_yflx_[xz]face, gr_zflx_[xy]face

#include "Simulation.h"

subroutine Grid_getFluxData_block(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, axis, isFluxDensity)
  use Grid_interface, ONLY : Grid_getCellFaceAreas

  use Grid_tile, ONLY : Grid_tile_t

  use Grid_data,      ONLY : gr_geometry

  use gr_specificData, ONLY : gr_xflx, gr_yflx, gr_zflx
  use gr_specificData, ONLY : gr_iloFl, gr_jloFl, gr_kloFl
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >=2
  use gr_specificData, ONLY : gr_xflx_yface, gr_yflx_xface
#if NDIM == 3
  use gr_specificData, ONLY : gr_xflx_zface, gr_yflx_zface, gr_zflx_xface, gr_zflx_yface
#endif
#endif
#endif
  use physicaldata, ONLY : flux_x, flux_y, flux_z, nfluxes
  use tree, ONLY : surr_blks

  implicit none

#include "constants.h"

  type(Grid_tile_t), intent(in) :: blockDesc
  integer,intent(in) :: lo(3)
  real,intent(OUT),dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxBufX, fluxBufY, fluxBufZ
  integer, intent(IN),optional :: axis
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

  real,pointer, dimension(:,:,:,:) :: fluxx,fluxy,fluxz

  integer :: blockID
  integer :: level              ! coarse level, for a block in which fluxes get corrected
  integer :: presVar
  integer :: sx,ex,sy,ey,sz,ez
  logical :: xtrue,ytrue,ztrue
  logical :: divideFluxx,divideFluxy,divideFluxz !whether to divide by area
  logical :: multFixed !whether to multiply by (1/2)**(NDIM-1)
#if NDIM == 1
  real,parameter :: scalInv = 1.0
#endif
#if NDIM == 2
  real,parameter :: scalInv = 0.5    ! 0.5**(NDIM-1)
#endif
#if NDIM == 3
  real,parameter :: scalInv = 0.25   ! 0.5**(NDIM-1)
#endif
  integer, dimension(MDIM) :: offs
  integer                  :: offx, offy, offz
  logical                  :: loCase2, hiCase2
  logical                  :: loCase4, hiCase4
  logical                  :: loCase24, hiCase24
  real,allocatable,target :: faceAreas(:,:,:)
  real,pointer            :: areaLeft(:,:,:)

  level = blockDesc % level
  blockID = blockDesc % id

  sx = NGUARD+1
  sy = NGUARD*K2D+1
  sz = NGUARD*K3D+1
#ifdef FIXEDBLOCKSIZE
  ex = NXB+NGUARD
  ey = NYB+NGUARD*K2D
  ez = NZB+NGUARD*K3D
#else
  ex = blockDesc % blkLimitsGC(HIGH,IAXIS)-NGUARD
  ey = blockDesc % blkLimitsGC(HIGH,JAXIS)-NGUARD*K2D
  ez = blockDesc % blkLimitsGC(HIGH,KAXIS)-NGUARD*K3D
#endif

  xtrue=.true.
  ytrue= (NDIM>1)
  ztrue= (NDIM>2)

  if(present(axis))then
     if (axis > 0 .AND. axis .NE. ALLDIR) then
        xtrue = (axis==IAXIS)
        ytrue = (axis==JAXIS)
        ztrue = (axis==KAXIS)
     end if
  end if

  select case (gr_geometry)
  case(CARTESIAN)
     divideFluxx = .false. ; divideFluxy = .false. ; divideFluxz = .false.
  case(SPHERICAL)
     divideFluxx = (NDIM>1); divideFluxy = .TRUE.  ; divideFluxz = .TRUE.
  case(POLAR)
     divideFluxx = .FALSE. ; divideFluxy = .FALSE. ; divideFluxz = .TRUE.
  case(CYLINDRICAL)
     divideFluxx = .FALSE. ; divideFluxy = .TRUE.  ; divideFluxz = .FALSE.
  end select

  multFixed = ((NDIM > 1) .AND. (gr_geometry .NE. CARTESIAN))

  if (nfluxes > 0) then
     fluxx(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxBufX ! fluxx,fluxy,fluxz use local (Paramesh) index counting
     fluxy(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxBufY ! fluxBuf[XYZ] use the global index convention (for the level)
     fluxz(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxBufZ
!!$     if (present(pressureSlots)) then
!!$        presP => pressureSlots
!!$     else
!!$        presP => presDefault
!!$     end if
     
     offs(:) = blockDesc%blkLimitsGC(LOW,1:MDIM) - 1
     offx = offs(IAXIS); offy = offs(JAXIS); offz = offs(KAXIS)

     if(xtrue) then
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >= 2
        fluxx(:,sx+2:ex-1,sy,sz:ez) =     gr_xflx_yface(:,3:NXB-1, 1     ,:,blockID)
        fluxx(:,sx+2:ex-1,ey,sz:ez) =     gr_xflx_yface(:,3:NXB-1, 2     ,:,blockID)
#if NDIM == 3
        fluxx(:,sx+2:ex-1,sy+1:ey-1,sz) = gr_xflx_zface(:,3:NXB-1,2:NYB-1,1,blockID)
        fluxx(:,sx+2:ex-1,sy+1:ey-1,ez) = gr_xflx_zface(:,3:NXB-1,2:NYB-1,2,blockID)
#endif
#endif
#endif

        ! With PARAMESH, there are four cases for the (same-level) neighbor in
        ! a certain direction of a *LEAF* block at refinement level blkLev:
        !
        ! Case | surr_blks(1,...) | surr_blks(3,...)  || Description: what's there?
        ! ======================= | ================= || ============================
        ! i.   !      <= -20      |  [ LEAF ? ]       || domain boundary (non-PERIODIC)
        ! ii.  |        -1        |  [  ignored  ]    || coarser block
        ! iii. | neighBlkID > 0   |  1  (LEAF)        || same refinement leaf block
        ! iv.  | neighBlkID > 0   |  2  (PARENT)      || finer blocks

        ! * We copy from flux_{x,y,z} in cases i,ii,iii,iv (i.e., all the above).
        ! * We divide by face areas if
        !   - other conditions are satisfied (geometry and direction), and
        !   - face area is nonzero, and
        !   - case is ii or iv.

        fluxx(:,sx,  sy:ey,sz:ez) = flux_x(:nfluxes,1,:,:,blockID)
        fluxx(:,ex+1,sy:ey,sz:ez) = flux_x(:nfluxes,2,:,:,blockID)
        fluxx(:,sx+1,sy:ey,sz:ez) = gr_xflx(:,1,:,:,blockID)
        fluxx(:,ex,  sy:ey,sz:ez) = gr_xflx(:,2,:,:,blockID)

        if (divideFluxx .OR. multFixed) then
           loCase2 = (surr_blks(1,1,1+K2D,1+K3D,blockID) == -1)
           hiCase2 = (surr_blks(1,3,1+K2D,1+K3D,blockID) == -1)
           loCase4 = (surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
                      surr_blks(3,1,1+K2D,1+K3D,blockID) == PARENT_BLK)
           hiCase4 = (surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
                      surr_blks(3,3,1+K2D,1+K3D,blockID) == PARENT_BLK)
           loCase24 = loCase2 .OR. loCase4
           hiCase24 = hiCase2 .OR. hiCase4

           if (divideFluxx .AND. (loCase24 .OR. hiCase24)) then
              allocate(faceAreas(offx+sx:offx+ex+1,offy+sy:offy+ey,offz+sz:offz+ez))
              areaLeft(sx:,sy:,sz:)  => faceAreas

              call Grid_getCellFaceAreas(IAXIS, level, &
                                         lbound(faceAreas), ubound(faceAreas), &
                                         faceAreas)
              do presVar = 1,nfluxes
                 if (loCase24) then
                    where (areaLeft(sx,sy:ey,sz:ez).NE.0.0) &
                         fluxx(presVar,sx,sy:ey,sz:ez) = fluxx(presVar,sx,sy:ey,sz:ez) / areaLeft(sx,sy:ey,sz:ez)
                 end if
                 if (hiCase24) then
                    fluxx(presVar,ex+1,sy:ey,sz:ez) = fluxx(presVar,ex+1,sy:ey,sz:ez) / areaLeft(ex+1,sy:ey,sz:ez)
                 end if
              end do
              deallocate(faceAreas)
           else if (loCase4 .OR. hiCase4) then
              do presVar = 1,nfluxes
                 if (loCase4) then
                    fluxx(presVar,sx,sy:ey,sz:ez) = fluxx(presVar,sx,sy:ey,sz:ez) * scalInv
                 end if
                 if (hiCase4) then
                    fluxx(presVar,ex+1,sy:ey,sz:ez) = fluxx(presVar,ex+1,sy:ey,sz:ez) * scalInv
                 end if
              end do
           end if
        end if
     end if
     
     if(ytrue) then
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >= 2
        fluxy(:,sx,sy+2:ey-1,sz:ez) = gr_yflx_xface(:,1,3:NYB-1,:,blockID)
        fluxy(:,ex,sy+2:ey-1,sz:ez) = gr_yflx_xface(:,2,3:NYB-1,:,blockID)
#if NDIM == 3
        fluxy(:,sx+1:ex-1,sy+2:ey-1,sz) = gr_yflx_zface(:,2:NXB-1,3:NYB-1,1,blockID)
        fluxy(:,sx+1:ex-1,sy+2:ey-1,ez) = gr_yflx_zface(:,2:NXB-1,3:NYB-1,2,blockID)
#endif
#endif
#endif

        fluxy(:,sx:ex,sy,  sz:ez) = flux_y(:nfluxes,:,1,:,blockID)
        fluxy(:,sx:ex,ey+1,sz:ez) = flux_y(:nfluxes,:,2,:,blockID)
        fluxy(:,sx:ex,sy+1,sz:ez) = gr_yflx(:,:,1,:,blockID)
        fluxy(:,sx:ex,ey,  sz:ez) = gr_yflx(:,:,2,:,blockID)
        
        
#if NDIM > 1
        if (divideFluxy .OR. multFixed) then
           loCase2 = (surr_blks(1,2,1,1+K3D,blockID) == -1)
           hiCase2 = (surr_blks(1,2,3,1+K3D,blockID) == -1)
           loCase4 = (surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
                      surr_blks(3,2,1,1+K3D,blockID) == PARENT_BLK)
           hiCase4 = (surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
                      surr_blks(3,2,3,1+K3D,blockID) == PARENT_BLK)
           loCase24 = loCase2 .OR. loCase4
           hiCase24 = hiCase2 .OR. hiCase4

           if (divideFluxy .AND. (loCase24 .OR. hiCase24)) then
              allocate(faceAreas(offx+sx:offx+ex,offy+sy:offy+ey+1,offz+sz:offz+ez))
              areaLeft(sx:,sy:,sz:)  => faceAreas

              call Grid_getCellFaceAreas(JAXIS, level, &
                   lbound(faceAreas), ubound(faceAreas), &
                   faceAreas)
              do presVar = 1,nfluxes
                 if (loCase24) then
                    where (areaLeft(sx:ex,sy,sz:ez).NE.0.0) &
                         fluxy(presVar,sx:ex,sy,sz:ez) = fluxy(presVar,sx:ex,sy,sz:ez) / areaLeft(sx:ex,sy,sz:ez)
                 end if
                 if (hiCase24) then
                    where (areaLeft(sx:ex,ey+1,sz:ez).NE.0.0) &
                         fluxy(presVar,sx:ex,ey+1,sz:ez) = fluxy(presVar,sx:ex,ey+1,sz:ez) / areaLeft(sx:ex,ey+1,sz:ez)
                 end if
              end do
              deallocate(faceAreas)
           else if (loCase4 .OR. hiCase4) then
              do presVar = 1,nfluxes
                 if (loCase4) then
                    fluxy(presVar,sx:ex,sy,sz:ez) = fluxy(presVar,sx:ex,sy,sz:ez) * scalInv
                 end if
                 if (hiCase24) then
                    fluxy(presVar,sx:ex,ey+1,sz:ez) = fluxy(presVar,sx:ex,ey+1,sz:ez) * scalInv
                 end if
              end do
           end if
        end if
#endif
     end if
     
     if(ztrue) then
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM == 3
        fluxz(:,sx,sy:ey,sz+2:ez-1) = gr_zflx_xface(:, 1,      : ,3:NZB-1,blockID)
        fluxz(:,ex,sy:ey,sz+2:ez-1) = gr_zflx_xface(:, 2,      : ,3:NZB-1,blockID)
        fluxz(:,sx+1:ex-1,sy,sz+2:ez-1) = gr_zflx_yface(:,2:NXB-1, 1 ,3:NZB-1,blockID)
        fluxz(:,sx+1:ex-1,ey,sz+2:ez-1) = gr_zflx_yface(:,2:NXB-1, 2 ,3:NZB-1,blockID)
#endif
#endif

        fluxz(:,sx:ex,sy:ey,sz  ) = flux_z(:nfluxes,:,:,1,blockID)
        fluxz(:,sx:ex,sy:ey,ez+1) = flux_z(:nfluxes,:,:,2,blockID)
        fluxz(:,sx:ex,sy:ey,sz+1) = gr_zflx(:,:,:,1,blockID)
        fluxz(:,sx:ex,sy:ey,ez  ) = gr_zflx(:,:,:,2,blockID)
#if NDIM > 2
        if (divideFluxz .OR. multFixed) then
           loCase2 = (surr_blks(1,2,2,1,blockID) == -1)
           hiCase2 = (surr_blks(1,2,2,3,blockID) == -1)
           loCase4 = (surr_blks(1,2,2,1,blockID) > 0 .AND. &
                      surr_blks(3,2,2,1,blockID) == PARENT_BLK)
           hiCase4 = (surr_blks(1,2,2,3,blockID) > 0 .AND. &
                      surr_blks(3,2,2,3,blockID) == PARENT_BLK)
           loCase24 = loCase2 .OR. loCase4
           hiCase24 = hiCase2 .OR. hiCase4

           if (divideFluxz .AND. (loCase24 .OR. hiCase24)) then
              allocate(faceAreas(offx+sx:offx+ex,offy+sy:offy+ey,offz+sz:offz+ez+1))
              areaLeft(sx:,sy:,sz:)  => faceAreas

              call Grid_getCellFaceAreas(KAXIS, level, &
                                         lbound(faceAreas), ubound(faceAreas), &
                                         faceAreas)
              do presVar = 1,nfluxes
                 !        where (areaLeft(sx:ex,sy:ey,sz).NE.0.0) &  ! should not happen in any supported geometry
                 if (loCase24) then
                    fluxz(presVar,sx:ex,sy:ey,sz) = fluxz(presVar,sx:ex,sy:ey,sz) / areaLeft(sx:ex,sy:ey,sz)
                 end if
                 if (hiCase24) then
                    fluxz(presVar,sx:ex,sy:ey,ez+1) = fluxz(presVar,sx:ex,sy:ey,ez+1) / areaLeft(sx:ex,sy:ey,ez+1)
                 end if
              end do
              deallocate(faceAreas)
           else if (loCase4 .OR. hiCase4) then
              do presVar = 1,nfluxes
                 if (loCase4) then
                    fluxz(presVar,sx:ex,sy:ey,sz) = fluxz(presVar,sx:ex,sy:ey,sz) * scalInv
                 end if
                 if (hiCase4) then
                    fluxz(presVar,sx:ex,sy:ey,ez+1) = fluxz(presVar,sx:ex,sy:ey,ez+1) * scalInv
                 end if
              end do
           end if
        end if
#endif
     end if
  end if

! At this point, border elements of fluxBufX have been updated through the fluxx pointer
! (unless an 'axis' argument has requested otherwise).
! Also, border elements of fluxBufY have been updated through the fluxy pointer
! (unless NDIM < 2 or an 'axis' argument has requested otherwise).
! Also, border elements of fluxBufZ have been updated through the fluxz pointer
! (unless NDIM < 3 or an 'axis' argument has requested otherwise).

end subroutine Grid_getFluxData_block
