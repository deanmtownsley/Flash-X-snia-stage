!!****if* source/Grid/GridMain/AMR/Paramesh4/Grid_correctFluxData_xtra
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
!!  Grid_correctFluxData_xtra
!!
!! SYNOPSIS
!!
!!  call Grid_correctFluxData_xtra(type(Grid_tile_t)(in) :: blockDesc,
!!                            real(IN)              :: scaleF,
!!                            real(INOUT),TARGET :: fluxBufX(size(fluxOldX,1),size(fluxOldX,2),size(fluxOldX,3),size(fluxOldX,4)),
!!                            real(INOUT),TARGET :: fluxBufY(size(fluxOldY,1),size(fluxOldY,2),size(fluxOldY,3),size(fluxOldY,4)),
!!                            real(INOUT),TARGET :: fluxBufZ(size(fluxOldZ,1),size(fluxOldZ,2),size(fluxOldZ,3),size(fluxOldZ,4)),
!!                            integer(in)           :: lo(MDIM),
!!                            real(IN)              :: scaleC,
!!                            real(in),   TARGET,CONTIGUOUS :: fluxOldX (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            real(in),   TARGET,CONTIGUOUS :: fluxOldY (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            real(in),   TARGET,CONTIGUOUS :: fluxOldZ (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            logical(IN), OPTIONAL  :: isFluxDensity)
!!
!!  call Grid_correctFluxData(type(Grid_tile_t)(in) :: blockDesc,
!!                            real(IN)              :: scaleF,
!!                            real(INOUT),TARGET :: fluxBufX(size(fluxOldX,1),size(fluxOldX,2),size(fluxOldX,3),size(fluxOldX,4)),
!!                            real(INOUT),TARGET :: fluxBufY(size(fluxOldY,1),size(fluxOldY,2),size(fluxOldY,3),size(fluxOldY,4)),
!!                            real(INOUT),TARGET :: fluxBufZ(size(fluxOldZ,1),size(fluxOldZ,2),size(fluxOldZ,3),size(fluxOldZ,4)),
!!                            integer(in)           :: lo(MDIM),
!!                            real(IN)              :: scaleC,
!!                            real(in),   TARGET,CONTIGUOUS :: fluxOldX (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            real(in),   TARGET,CONTIGUOUS :: fluxOldY (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            real(in),   TARGET,CONTIGUOUS :: fluxOldZ (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            logical(IN), OPTIONAL  :: isFluxDensity)
!!
!! DESCRIPTION
!!
!!   Correct data in flux arrays by replacing fluxes in certain locations with
!!   a linear combination of data from a higher refinement level and
!!   a previously computed approximation of coarse fluxes.
!!
!!     fluxBuf  :=  scaleF*"communicated fine fluxes" + scaleC*fluxOld  AT coarse side of f/c bdry;
!!              unmodified                                                  ELSEWHERE.
!!
!!   By proper choice of the sign of the scaling factors, the previously
!!   computed coarse fluxes can effectively be subtracted from the
!!   fluxes that come from a higher refinement level.
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
!!   scaleF   : coefficient for fluxes computed at finer resolution
!!
!!   fluxBufX : fluxes for IAXIS direction
!!
!!   fluxBufY : fluxes for JAXIS direction
!!
!!   fluxBufZ : fluxes for KAXIS direction
!!
!!   lo :   lower bounds for the spatial indices of the flux buffers
!!
!!   scaleC   : coefficient for previously computed coarse fluxes
!!
!!   fluxOldX : previously computed coarse fluxes for IAXIS direction
!!
!!   fluxOldY : previously computed coarse fluxes for JAXIS direction
!!
!!   fluxOldZ : previously computed coarse fluxes for KAXIS direction
!!
!!   isFluxDensity : are the fluxes actually fluxes or flux densities?
!!
!! NOTES
!!
!!   This subroutine is available under the generic name Grid_correctFluxData
!!   as well as under the specifc name Grid_correctFluxData_xtra.
!!   The calling code should refer to the appropriate name in a statement
!!      use Grid_interface, ONLY: ...
!!
!!   The arrays fluxOldX, fluxOldY, fluxOldZ are subject to index reordering.
!!   The arrays fluxBufX, fluxBufY, fluxBufZ are index-reorder indirectly,
!!   by having bounds that depend  on the shape on the corresponding
!!   fluxOldX, fluxOldY, fluxOldZ arrays.
!!
!!   flux buffer arrays and flux correction arrays should contain
!!   space for fluxes of all valid cells in the block, excluding guard
!!   cells.
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

!!REORDER(4): fluxOld[XYZ]
!!REORDER(5): flux_[xyz]
!!REORDER(4): fluxx,fluxy,fluxz,prevCoarse[XYZ]

#include "Simulation.h"

subroutine Grid_correctFluxData_xtra(blockDesc, scaleF, fluxBufX,fluxBufY,fluxBufZ, lo, &
                                                scaleC, fluxOldX,fluxOldY,fluxOldZ,     &
                                                isFluxDensity)
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
  real,intent(in)    :: scaleF,scaleC
  real,intent(in)   ,dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxOldX, fluxOldY, fluxOldZ
  CONTIGUOUS_FSTMT(fluxOldX)
  CONTIGUOUS_FSTMT(fluxOldY)
  CONTIGUOUS_FSTMT(fluxOldZ)
  real,INTENT(INOUT),dimension(size(fluxOldX,1),size(fluxOldX,2),size(fluxOldX,3),size(fluxOldX,4)),TARGET :: fluxBufX
  real,INTENT(INOUT),dimension(size(fluxOldY,1),size(fluxOldY,2),size(fluxOldY,3),size(fluxOldY,4)),TARGET :: fluxBufY
  real,INTENT(INOUT),dimension(size(fluxOldZ,1),size(fluxOldZ,2),size(fluxOldZ,3),size(fluxOldZ,4)),TARGET :: fluxBufZ
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

  real,pointer, dimension(:,:,:,:) :: fluxx,fluxy,fluxz
  real,pointer, dimension(:,:,:,:) :: prevCoarseX,prevCoarseY,prevCoarseZ
  CONTIGUOUS_FSTMT(fluxx)
  CONTIGUOUS_FSTMT(fluxy)
  CONTIGUOUS_FSTMT(fluxz)
  CONTIGUOUS_FSTMT(prevCoarseX)
  CONTIGUOUS_FSTMT(prevCoarseY)
  CONTIGUOUS_FSTMT(prevCoarseZ)

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
     fluxy(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxBufY ! Actual arguments for fluxBuf[XYZ] effectively use the same index convention
     fluxz(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxBufZ ! as fluxOld[XYZ].

     prevCoarseX(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxOldX ! prevCoarseX,prevCoarseY,prevCoarseZ use local (Paramesh) index counting
     prevCoarseY(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxOldY ! fluxOld[XYZ] use the global index convention (for the level)
     prevCoarseZ(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxOldZ
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
        ! * We divide the linear combination of flux_{x,y,z} and tflux_{x,y,z}
        !   by face areas if
        !   - other conditions are satisfied (geometry and direction), and
        !   - face area is nonzero, and
        !   - case is iv.

        loCase4 = (surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,1,1+K2D,1+K3D,blockID) == PARENT_BLK)
        hiCase4 = (surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,3,1+K2D,1+K3D,blockID) == PARENT_BLK)

        if (loCase4) fluxx(:,sx,  sy:ey,sz:ez) = scaleF*flux_x(:nfluxes,1,:,:,blockID)
        if (hiCase4) fluxx(:,ex+1,sy:ey,sz:ez) = scaleF*flux_x(:nfluxes,2,:,:,blockID)

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
           if (loCase4) fluxx(:,sx,  sy:ey,sz:ez) = fluxx(:,sx,  sy:ey,sz:ez) + scaleC*prevCoarseX(:,sx,  sy:ey,sz:ez)
           if (hiCase4) fluxx(:,ex+1,sy:ey,sz:ez) = fluxx(:,ex+1,sy:ey,sz:ez) + scaleC*prevCoarseX(:,ex+1,sy:ey,sz:ez)
        end if
     end if

     if(ytrue) then

#if NDIM > 1
        loCase4 = (surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,2,1,1+K3D,blockID) == PARENT_BLK)
        hiCase4 = (surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
                   surr_blks(3,2,3,1+K3D,blockID) == PARENT_BLK)
        if (loCase4) fluxy(:,sx:ex,sy,  sz:ez) = scaleF*flux_y(:nfluxes,:,1,:,blockID)
        if (hiCase4) fluxy(:,sx:ex,ey+1,sz:ez) = scaleF*flux_y(:nfluxes,:,2,:,blockID)


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
           if (loCase4) fluxy(:,sx:ex,sy,  sz:ez) = fluxy(:,sx:ex,sy,  sz:ez) + scaleC*prevCoarseY(:,sx:ex,sy,  sz:ez)
           if (hiCase4) fluxy(:,sx:ex,ey+1,sz:ez) = fluxy(:,sx:ex,ey+1,sz:ez) + scaleC*prevCoarseY(:,sx:ex,ey+1,sz:ez)
        end if
#endif
     end if

#if NDIM > 2
     if(ztrue) then

        loCase4 = (surr_blks(1,2,2,1,blockID) > 0 .AND. &
                   surr_blks(3,2,2,1,blockID) == PARENT_BLK)
        hiCase4 = (surr_blks(1,2,2,3,blockID) > 0 .AND. &
                   surr_blks(3,2,2,3,blockID) == PARENT_BLK)
        if (loCase4) fluxz(:,sx:ex,sy:ey,sz  ) = scaleF*flux_z(:nfluxes,:,:,1,blockID)
        if (hiCase4) fluxz(:,sx:ex,sy:ey,ez+1) = scaleF*flux_z(:nfluxes,:,:,2,blockID)

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
           if (loCase4) fluxz(:,sx:ex,sy:ey,sz  ) = fluxz(:,sx:ex,sy:ey,sz  ) + scaleC*prevCoarseZ(:,sx:ex,sy:ey,sz  )
           if (hiCase4) fluxz(:,sx:ex,sy:ey,ez+1) = fluxz(:,sx:ex,sy:ey,ez+1) + scaleC*prevCoarseZ(:,sx:ex,sy:ey,ez+1)
        end if
     end if
#endif
  end if

end subroutine Grid_correctFluxData_xtra
