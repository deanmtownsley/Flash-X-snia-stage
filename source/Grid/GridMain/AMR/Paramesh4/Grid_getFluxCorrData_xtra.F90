!!****if* source/Grid/GridMain/AMR/Paramesh4/Grid_getFluxCorrData_xtra
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
!!  Grid_getFluxCorrData_xtra
!!
!! SYNOPSIS
!!
!!  call Grid_getFluxCorrData_xtra(type(Grid_tile_t)(in) :: blockDesc,
!!                                 real(IN), TARGET,CONTIGUOUS :: fluxBufX(:,lo(1): ,lo(2): ,lo(3): ),
!!                                 real(IN), TARGET,CONTIGUOUS :: fluxBufY(:,lo(1): ,lo(2): ,lo(3): ),
!!                                 real(IN), TARGET,CONTIGUOUS :: fluxBufZ(:,lo(1): ,lo(2): ,lo(3): ),
!!                                 integer(in)           :: lo(3),
!!                                 real(OUT),TARGET,CONTIGUOUS :: fluxCorrX(:,lo(1): ,lo(2): ,lo(3): ),
!!                                 real(OUT),TARGET,CONTIGUOUS :: fluxCorrY(:,lo(1): ,lo(2): ,lo(3): ),
!!                                 real(OUT),TARGET,CONTIGUOUS :: fluxCorrZ(:,lo(1): ,lo(2): ,lo(3): ),
!!                                 logical(IN), OPTIONAL :: isFluxDensity(:))
!!
!! DESCRIPTION
!!
!!   Get flux corrections from semipermanent flux storage (SPFS).
!!
!!     fluxCorr :=  "communicated fine fluxes" - fluxBuf  AT coarse side of f/c bdry;
!!              :=  0.0                                   ELSEWHERE.
!!
!!   Flux corrections are returned in fluxCorr[XYZ] arguments.
!!   The arguments fluxBuf[XYZ] are input only and represent coarse
!!   fluxes, i.e., the ones for which corrections need to be computed.
!!
!!   Only fluxes at locations that represent the coarse side of fine/coarse
!!   block boundaries can hold nonzero flux correction data on return.
!!   Other elements of the flux buffers are set to zero if they represent
!!   faces of any cells that touch a block boundary; data faces of cells
!!   that are farther away from a block boundary are left undefined.
!!
!! ARGUMENTS
!!
!!   blockDesc : describes the current block.
!!               Note that this should be a full block, not a tile representing
!!               a partial block.
!!
!!   fluxBufX :  buffer for fluxes in IAXIS-direction
!!
!!   fluxBufY :  buffer for fluxes in JAXIS-direction; ignored if NDIM < 2
!!
!!   fluxBufZ :  buffer for fluxes in KAXIS-direction; ignored if NDIM < 3
!!
!!   lo :        lower bounds for the spatial indices of the flux buffers
!!
!!   fluxCorrX : flux correction (difference) for IAXIS direction
!!
!!   fluxCorrY : flux correction (difference) for JAXIS direction;
!!               left undefined if NDIM < 2.
!!
!!   fluxCorrZ : flux correction (difference) for KAXIS direction;
!!               left undefined if NDIM < 3.
!!
!!   isFluxDensity : indicates, for each flux component, whether the component
!!                   is a flux proper (if TRUE) or a flux density (otherwise).
!!                   This may be either removed, or changed into a scalar flag,
!!                   later.
!!
!! NOTES
!!
!!   The arrays fluxBufX, fluxBufY, fluxBufZ are subject to index reordering.
!!   The arrays fluxCorrX, fluxCorrY, fluxCorrZ are subject to index reordering.
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

!!REORDER(4): fluxBuf[XYZ],fluxCorr[XYZ]
!!REORDER(5): flux_[xyz]
!!REORDER(4): fluxx,fluxy,fluxz,corrx,corry,corrz

#include "Simulation.h"

subroutine Grid_getFluxCorrData_xtra(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, fluxCorrX,fluxCorrY,fluxCorrZ, isFluxDensity)
  use Grid_interface, ONLY : Grid_getCellFaceAreas

  use Grid_tile, ONLY : Grid_tile_t

  use Grid_data,      ONLY : gr_geometry

  use gr_specificData, ONLY : gr_iloFl, gr_jloFl, gr_kloFl
  use physicaldata, ONLY : flux_x, flux_y, flux_z, nfluxes
  use tree, ONLY : surr_blks

  implicit none

#include "constants.h"
#include "FortranLangFeatures.fh"

  type(Grid_tile_t), intent(in) :: blockDesc
  integer,intent(in) :: lo(3)
  real,intent(in)   ,dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxBufX, fluxBufY, fluxBufZ
  real,intent(OUT)  ,dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxCorrX, fluxCorrY, fluxCorrZ
  CONTIGUOUS_FSTMT(fluxBufX)
  CONTIGUOUS_FSTMT(fluxBufY)
  CONTIGUOUS_FSTMT(fluxBufZ)
  CONTIGUOUS_FSTMT(fluxCorrX)
  CONTIGUOUS_FSTMT(fluxCorrY)
  CONTIGUOUS_FSTMT(fluxCorrZ)
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

  real,pointer, dimension(:,:,:,:) :: fluxx,fluxy,fluxz
  real,pointer, dimension(:,:,:,:) :: corrx,corry,corrz
  CONTIGUOUS_FSTMT(fluxx)
  CONTIGUOUS_FSTMT(fluxy)
  CONTIGUOUS_FSTMT(fluxz)
  CONTIGUOUS_FSTMT(corrx)
  CONTIGUOUS_FSTMT(corry)
  CONTIGUOUS_FSTMT(corrz)

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
  logical                  :: loCase4, hiCase4
  real,allocatable,target :: faceAreas(:,:,:)
  real,pointer            :: areaLeft(:,:,:)
  CONTIGUOUS_FSTMT(areaLeft)

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

     corrx(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxCorrX ! corrx,corry,corrz use local (Paramesh) index counting
     corry(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxCorrY ! fluxCorr[XYZ] use the global index convention (for the level)
     corrz(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxCorrZ
!!$     if (present(pressureSlots)) then
!!$        presP => pressureSlots
!!$     else
!!$        presP => presDefault
!!$     end if

     offs(:) = blockDesc%blkLimitsGC(LOW,1:MDIM) - 1
     offx = offs(IAXIS); offy = offs(JAXIS); offz = offs(KAXIS)

     if(xtrue) then

#if NDIM >= 2
        corrx(:,sx+2:ex-1,sy,sz:ez) = 0.0
        corrx(:,sx+2:ex-1,ey,sz:ez) = 0.0
#if NDIM == 3
        corrx(:,sx+2:ex-1,sy+1:ey-1,sz) = 0.0
        corrx(:,sx+2:ex-1,sy+1:ey-1,ez) = 0.0
#endif
#endif

        corrx(:,sx:sx+1,sy:ey,sz:ez) = 0.0
        corrx(:,ex:ex+1,sy:ey,sz:ez) = 0.0

        ! With PARAMESH, there are four cases for the (same-level) neighbor in
        ! a certain direction of a *LEAF* block at refinement level blkLev:
        !
        ! Case | surr_blks(1,...) | surr_blks(3,...)  || Description: what's there?
        ! ======================= | ================= || ============================
        ! i.   !      <= -20      |  [ LEAF ? ]       || domain boundary (non-PERIODIC)
        ! ii.  |        -1        |  [  ignored  ]    || coarser block
        ! iii. | neighBlkID > 0   |  1  (LEAF)        || same refinement leaf block
        ! iv.  | neighBlkID > 0   |  2  (PARENT)      || finer blocks

        ! * We copy from flux_{x,y,z} in case iv.
        ! * We divide the difference between flux_{x,y,z} and fluxBuf{X,Y,Z}
        !   by face areas if
        !   - other conditions are satisfied (geometry and direction), and
        !   - face area is nonzero, and
        !   - case is iv.

        loCase4 = (surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,1,1+K2D,1+K3D,blockID) == PARENT_BLK)
        hiCase4 = (surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,3,1+K2D,1+K3D,blockID) == PARENT_BLK)

        if (loCase4) corrx(:,sx,  sy:ey,sz:ez) = flux_x(:nfluxes,1,:,:,blockID)
        if (hiCase4) corrx(:,ex+1,sy:ey,sz:ez) = flux_x(:nfluxes,2,:,:,blockID)

        if (loCase4 .OR. hiCase4) then

           if (divideFluxx) then
              allocate(faceAreas(offx+sx:offx+ex+1,offy+sy:offy+ey,offz+sz:offz+ez))
              areaLeft(sx:,sy:,sz:)  => faceAreas

              call Grid_getCellFaceAreas(IAXIS, level, &
                                         lbound(faceAreas), ubound(faceAreas), &
                                         faceAreas)
              do presVar = 1,nfluxes
                 if (loCase4) then
                    where (areaLeft(sx,sy:ey,sz:ez).NE.0.0) &
                         corrx(presVar,sx,sy:ey,sz:ez) = corrx(presVar,sx,sy:ey,sz:ez) / areaLeft(sx,sy:ey,sz:ez)
                 end if
                 if (hiCase4) then
                    corrx(presVar,ex+1,sy:ey,sz:ez) = corrx(presVar,ex+1,sy:ey,sz:ez) / areaLeft(ex+1,sy:ey,sz:ez)
                 end if
              end do
              deallocate(faceAreas)
           else if (multFixed) then
              do presVar = 1,nfluxes
                 if (loCase4) then
                    corrx(presVar,sx,sy:ey,sz:ez) = corrx(presVar,sx,sy:ey,sz:ez) * scalInv
                 end if
                 if (hiCase4) then
                    corrx(presVar,ex+1,sy:ey,sz:ez) = corrx(presVar,ex+1,sy:ey,sz:ez) * scalInv
                 end if
              end do
           end if
           if (loCase4) corrx(:,sx,  sy:ey,sz:ez) = corrx(:,sx,  sy:ey,sz:ez) - fluxx(:,sx,  sy:ey,sz:ez)
           if (hiCase4) corrx(:,ex+1,sy:ey,sz:ez) = corrx(:,ex+1,sy:ey,sz:ez) - fluxx(:,ex+1,sy:ey,sz:ez)
        end if
     end if

     if(ytrue) then

#if NDIM >= 2
        corry(:,sx,sy+2:ey-1,sz:ez) = 0.0
        corry(:,ex,sy+2:ey-1,sz:ez) = 0.0
#if NDIM == 3
        corry(:,sx+1:ex-1,sy+2:ey-1,sz) = 0.0
        corry(:,sx+1:ex-1,sy+2:ey-1,ez) = 0.0
#endif
#endif

        corry(:,sx:ex,sy:sy+1,sz:ez) = 0.0
        corry(:,sx:ex,ey:ey+1,sz:ez) = 0.0

#if NDIM > 1
        loCase4 = (surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,2,1,1+K3D,blockID) == PARENT_BLK)
        hiCase4 = (surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,2,3,1+K3D,blockID) == PARENT_BLK)
        if (loCase4) corry(:,sx:ex,sy,  sz:ez) = flux_y(:nfluxes,:,1,:,blockID)
        if (hiCase4) corry(:,sx:ex,ey+1,sz:ez) = flux_y(:nfluxes,:,2,:,blockID)


        if (loCase4 .OR. hiCase4) then

           if (divideFluxy) then
              allocate(faceAreas(offx+sx:offx+ex,offy+sy:offy+ey+1,offz+sz:offz+ez))
              areaLeft(sx:,sy:,sz:)  => faceAreas

              call Grid_getCellFaceAreas(JAXIS, level, &
                   lbound(faceAreas), ubound(faceAreas), &
                   faceAreas)
              do presVar = 1,nfluxes
                 if (loCase4) then
                    where (areaLeft(sx:ex,sy,sz:ez).NE.0.0) &
                         corry(presVar,sx:ex,sy,sz:ez) = corry(presVar,sx:ex,sy,sz:ez) / areaLeft(sx:ex,sy,sz:ez)
                 end if
                 if (hiCase4) then
                    where (areaLeft(sx:ex,ey+1,sz:ez).NE.0.0) &
                         corry(presVar,sx:ex,ey+1,sz:ez) = corry(presVar,sx:ex,ey+1,sz:ez) / areaLeft(sx:ex,ey+1,sz:ez)
                 end if
              end do
              deallocate(faceAreas)
           else if (multFixed) then
              do presVar = 1,nfluxes
                 if (loCase4) then
                    corry(presVar,sx:ex,sy,sz:ez) = corry(presVar,sx:ex,sy,sz:ez) * scalInv
                 end if
                 if (hiCase4) then
                    corry(presVar,sx:ex,ey+1,sz:ez) = corry(presVar,sx:ex,ey+1,sz:ez) * scalInv
                 end if
              end do
           end if
           if (loCase4) corry(:,sx:ex,sy,  sz:ez) = corry(:,sx:ex,sy,  sz:ez) - fluxy(:,sx:ex,sy,  sz:ez)
           if (hiCase4) corry(:,sx:ex,ey+1,sz:ez) = corry(:,sx:ex,ey+1,sz:ez) - fluxy(:,sx:ex,ey+1,sz:ez)
        end if
#endif
     end if

     if(ztrue) then

#if NDIM > 2
        corrz(:,sx       ,sy:ey,sz+2:ez-1) = 0.0
        corrz(:,ex       ,sy:ey,sz+2:ez-1) = 0.0
        corrz(:,sx+1:ex-1,sy   ,sz+2:ez-1) = 0.0
        corrz(:,sx+1:ex-1,ey   ,sz+2:ez-1) = 0.0

        corrz(:,sx:ex    ,sy:ey,  sz:sz+1) = 0.0
        corrz(:,sx:ex    ,sy:ey,  ez:ez+1) = 0.0

        loCase4 = (surr_blks(1,2,2,1,blockID) > 0 .AND. &
                   surr_blks(3,2,2,1,blockID) == PARENT_BLK)
        hiCase4 = (surr_blks(1,2,2,3,blockID) > 0 .AND. &
                   surr_blks(3,2,2,3,blockID) == PARENT_BLK)
        if (loCase4) corrz(:,sx:ex,sy:ey,sz  ) = flux_z(:nfluxes,:,:,1,blockID)
        if (hiCase4) corrz(:,sx:ex,sy:ey,ez+1) = flux_z(:nfluxes,:,:,2,blockID)

        if (loCase4 .OR. hiCase4) then

           if (divideFluxz) then
              allocate(faceAreas(offx+sx:offx+ex,offy+sy:offy+ey,offz+sz:offz+ez+1))
              areaLeft(sx:,sy:,sz:)  => faceAreas

              call Grid_getCellFaceAreas(KAXIS, level, &
                                         lbound(faceAreas), ubound(faceAreas), &
                                         faceAreas)
              do presVar = 1,nfluxes
                 !        where (areaLeft(sx:ex,sy:ey,sz).NE.0.0) &  ! should not happen in any supported geometry
                 if (loCase4) then
                    corrz(presVar,sx:ex,sy:ey,sz) = corrz(presVar,sx:ex,sy:ey,sz) / areaLeft(sx:ex,sy:ey,sz)
                 end if
                 if (hiCase4) then
                    corrz(presVar,sx:ex,sy:ey,ez+1) = corrz(presVar,sx:ex,sy:ey,ez+1) / areaLeft(sx:ex,sy:ey,ez+1)
                 end if
              end do
              deallocate(faceAreas)
           else if (multFixed) then
              do presVar = 1,nfluxes
                 if (loCase4) then
                    corrz(presVar,sx:ex,sy:ey,sz) = corrz(presVar,sx:ex,sy:ey,sz) * scalInv
                 end if
                 if (hiCase4) then
                    corrz(presVar,sx:ex,sy:ey,ez+1) = corrz(presVar,sx:ex,sy:ey,ez+1) * scalInv
                 end if
              end do
           end if
           if (loCase4) corrz(:,sx:ex,sy:ey,sz  ) = corrz(:,sx:ex,sy:ey,sz  ) - fluxz(:,sx:ex,sy:ey,sz  )
           if (hiCase4) corrz(:,sx:ex,sy:ey,ez+1) = corrz(:,sx:ex,sy:ey,ez+1) - fluxz(:,sx:ex,sy:ey,ez+1)
        end if
#endif
     end if
  end if

end subroutine Grid_getFluxCorrData_xtra
