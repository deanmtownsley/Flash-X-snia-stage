!!****if* source/Grid/GridMain/AMR/Paramesh4/Grid_correctFluxData
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
!!  Grid_correctFluxData
!!
!! SYNOPSIS
!!
!!  call Grid_correctFluxData(type(Grid_tile_t)(in) :: blockDesc,
!!                            real(INOUT),TARGET,CONTIGUOUS :: fluxBufX (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            real(INOUT),TARGET,CONTIGUOUS :: fluxBufY (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            real(INOUT),TARGET,CONTIGUOUS :: fluxBufZ (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            integer(in)           :: lo(MDIM),
!!                            logical(IN), OPTIONAL :: isFluxDensity)
!!
!! DESCRIPTION
!!
!!   Correct data in flux arrays by replacing fluxes in certain locations with
!!   data from a higher refinement level.
!!
!!     fluxBuf  :=  "communicated fine fluxes"  AT        coarse side of f/c bdry;
!!              unmodified                      ELSEWHERE.
!!
!!   Finer-level data (where needed) must have been stored to SPFS,
!!   typically by calling Grid_putFluxData_block on relevant
!!   neigboring blocks, and communication must have been triggered,
!!   typically by calling Grid_communicateFluxes, before this
!!   interface is invoked for a block.
!!
!!   Only fluxes at locations that represent the coarse side of fine/coarse
!!   block boundaries are modified, other elements of the flux buffers are
!!   left unmodified by calling this interface.
!!
!! ARGUMENTS
!!
!!   blockDesc : descriptor for one block. !!DEV: can it be a proper tile?
!!
!!   fluxBufX : fluxes for IAXIS direction
!!
!!   fluxBufY : fluxes for JAXIS direction
!!
!!   fluxBufZ : fluxes for KAXIS direction
!!
!!   lo :   lower bounds for the spatial indices of the flux buffers
!!
!!   isFluxDensity : are the fluxes actually fluxes or flux densities?
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
!!   Hydro
!!
!!***

!!REORDER(4): fluxBuf[XYZ]
!!REORDER(5): flux_[xyz]
!!REORDER(4): fluxx,fluxy,fluxz

#include "Simulation.h"

subroutine Grid_correctFluxData(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, isFluxDensity)
  use Grid_interface, ONLY : Grid_getCellFaceAreas

  use Grid_tile, ONLY : Grid_tile_t

  use Grid_data,      ONLY : gr_geometry

  use gr_specificData, ONLY : gr_iloFl, gr_jloFl, gr_kloFl
  use physicaldata, ONLY : flux_x, flux_y, flux_z, nfluxes
  use tree, ONLY : surr_blks

  implicit none

#include "FortranLangFeatures.fh"
#include "constants.h"

  type(Grid_tile_t), intent(in) :: blockDesc
  integer,intent(in) :: lo(3)
  real,intent(INOUT),dimension(: , lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxBufX, fluxBufY, fluxBufZ
  CONTIGUOUS_FSTMT(fluxBufX)
  CONTIGUOUS_FSTMT(fluxBufY)
  CONTIGUOUS_FSTMT(fluxBufZ)
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

  real,pointer, dimension(:,:,:,:) :: fluxx,fluxy,fluxz
  CONTIGUOUS_FSTMT(fluxx)
  CONTIGUOUS_FSTMT(fluxy)
  CONTIGUOUS_FSTMT(fluxz)

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
!!$     if (present(pressureSlots)) then
!!$        presP => pressureSlots
!!$     else
!!$        presP => presDefault
!!$     end if

     offs(:) = blockDesc%blkLimitsGC(LOW,1:MDIM) - 1
     offx = offs(IAXIS); offy = offs(JAXIS); offz = offs(KAXIS)

     if(xtrue) then

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
        ! * We divide by face areas if
        !   - other conditions are satisfied (geometry and direction), and
        !   - face area is nonzero, and
        !   - case is iv.

        loCase4 = (surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,1,1+K2D,1+K3D,blockID) == PARENT_BLK)
        hiCase4 = (surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,3,1+K2D,1+K3D,blockID) == PARENT_BLK)

        if (loCase4) fluxx(:,sx,  sy:ey,sz:ez) = flux_x(:nfluxes,1,:,:,blockID)
        if (hiCase4) fluxx(:,ex+1,sy:ey,sz:ez) = flux_x(:nfluxes,2,:,:,blockID)

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
                         fluxx(presVar,sx,sy:ey,sz:ez) = fluxx(presVar,sx,sy:ey,sz:ez) / areaLeft(sx,sy:ey,sz:ez)
                 end if
                 if (hiCase4) then
                    fluxx(presVar,ex+1,sy:ey,sz:ez) = fluxx(presVar,ex+1,sy:ey,sz:ez) / areaLeft(ex+1,sy:ey,sz:ez)
                 end if
              end do
              deallocate(faceAreas)
           else if (multFixed) then
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

#if NDIM > 1
        loCase4 = (surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,2,1,1+K3D,blockID) == PARENT_BLK)
        hiCase4 = (surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,2,3,1+K3D,blockID) == PARENT_BLK)
        if (loCase4) fluxy(:,sx:ex,sy,  sz:ez) = flux_y(:nfluxes,:,1,:,blockID)
        if (hiCase4) fluxy(:,sx:ex,ey+1,sz:ez) = flux_y(:nfluxes,:,2,:,blockID)


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
                         fluxy(presVar,sx:ex,sy,sz:ez) = fluxy(presVar,sx:ex,sy,sz:ez) / areaLeft(sx:ex,sy,sz:ez)
                 end if
                 if (hiCase4) then
                    where (areaLeft(sx:ex,ey+1,sz:ez).NE.0.0) &
                         fluxy(presVar,sx:ex,ey+1,sz:ez) = fluxy(presVar,sx:ex,ey+1,sz:ez) / areaLeft(sx:ex,ey+1,sz:ez)
                 end if
              end do
              deallocate(faceAreas)
           else if (multFixed) then
              do presVar = 1,nfluxes
                 if (loCase4) then
                    fluxy(presVar,sx:ex,sy,sz:ez) = fluxy(presVar,sx:ex,sy,sz:ez) * scalInv
                 end if
                 if (hiCase4) then
                    fluxy(presVar,sx:ex,ey+1,sz:ez) = fluxy(presVar,sx:ex,ey+1,sz:ez) * scalInv
                 end if
              end do
           end if
        end if
#endif
     end if

#if NDIM > 2
     if(ztrue) then

        loCase4 = (surr_blks(1,2,2,1,blockID) > 0 .AND. &
                   surr_blks(3,2,2,1,blockID) == PARENT_BLK)
        hiCase4 = (surr_blks(1,2,2,3,blockID) > 0 .AND. &
                   surr_blks(3,2,2,3,blockID) == PARENT_BLK)
        if (loCase4) fluxz(:,sx:ex,sy:ey,sz  ) = flux_z(:nfluxes,:,:,1,blockID)
        if (hiCase4) fluxz(:,sx:ex,sy:ey,ez+1) = flux_z(:nfluxes,:,:,2,blockID)

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
                    fluxz(presVar,sx:ex,sy:ey,sz) = fluxz(presVar,sx:ex,sy:ey,sz) / areaLeft(sx:ex,sy:ey,sz)
                 end if
                 if (hiCase4) then
                    fluxz(presVar,sx:ex,sy:ey,ez+1) = fluxz(presVar,sx:ex,sy:ey,ez+1) / areaLeft(sx:ex,sy:ey,ez+1)
                 end if
              end do
              deallocate(faceAreas)
           else if (multFixed) then
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
     end if
#endif
  end if

end subroutine Grid_correctFluxData
