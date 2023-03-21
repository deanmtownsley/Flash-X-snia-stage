!!****if* source/Grid/GridMain/AMR/Paramesh4/Grid_putFluxData_block
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  Grid_putFluxData_block
!!
!! SYNOPSIS
!!
!!  call Grid_putFluxData_block(type(Grid_tile_t)(in) :: blockDesc,
!!                              real(in),Contiguous,TARGET :: fluxBufX(:, lo(1): ,lo(2): ,lo(3): ),
!!                              real(in),Contiguous,TARGET :: fluxBufY(:, lo(1): ,lo(2): ,lo(3): ),
!!                              real(in),Contiguous,TARGET :: fluxBufZ(:, lo(1): ,lo(2): ,lo(3): ),
!!                              integer(in)           :: lo(3),
!!                              logical(IN), OPTIONAL :: add,
!!                              logical(IN), OPTIONAL :: isFluxDensity)
!!
!! DESCRIPTION
!!
!!   Save fluxes passed in as arguments in semipermanent flux storage (SPFS).
!!
!!   For the PARAMESH Grid implementation, SPFS translates to the private arrays
!!   flux_x, flux_y, flux_z, and maybe additional arrays like gr_[xyz]flx,
!!   gr_xflx_[yz]face, gr_yflx_[xz]face, gr_zflx_[xy]face.
!!
!!   This implementation uses the traditional PARAMESH 4 storage for SPFS.
!!
!! ARGUMENTS
!!
!!   blockDesc : describes the current block.
!!               Note that this should be a full block, not a tile representing
!!               a partial block.
!!
!!   fluxBufX :  fluxes in IAXIS-direction
!!
!!   fluxBufY :  fluxes in JAXIS-direction; ignored if NDIM < 2
!!
!!   fluxBufZ :  fluxes in KAXIS-direction; ignored if NDIM < 3
!!
!!   lo :        lower bounds for the spatial indices of the flux buffers
!!
!!   add :       whether to add or override.
!!               If this argument is present in a call but isFluxDensity is not,
!!               INDEPENDENT OF WHETHER THE VALUE IS .FALSE. OR .TRUE.,
!!               it is assumed that in the non-Cartesian geometry the "fluxes"
!!               in the fluxBuf? arguments have already been multiplied by
!!               area factors; that is, they really are fluxes and not flux
!!               densities.
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
!!   This interface does not require level-wide fluxes to be allocated.
!!
!!   SPFS means semi-permanent flux storage. When using a Grid
!!   implementation based on PARAMESH, SPFS is provided by arrays
!!   flux_x, flux_y, flux_z private to PARAMESH in conjunction with
!!   additional arrays like gr_[xyz]flx and gr_xflx_[yz]face,gr_yflx_[xz]face,
!!   gr_zflx_[xy]face that are private to the Paramesh4 immplementation of
!!   the Grid.
!!
!! SEE ALSO
!!
!!   Grid_communicateFluxes
!!   Grid_correctFluxData
!!   Hydro
!!***

!!REORDER(4): fluxBuf[XYZ]
!!REORDER(5): flux_[xyz],gr_[xyz]flx
!!REORDER(4): fluxx,fluxy,fluxz
!!REORDER(5): gr_xflx_[yz]face, gr_yflx_[xz]face, gr_zflx_[xy]face

#include "Simulation.h"
subroutine Grid_putFluxData_block(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, add, isFluxDensity)
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
  use tree, ONLY : surr_blks, nodetype, lrefine_max

  implicit none

#include "constants.h"
#include "FortranLangFeatures.fh"

  type(Grid_tile_t), intent(in) :: blockDesc
  integer,intent(in) :: lo(3)
  real,intent(in),dimension(:, lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxBufX,fluxBufY,fluxBufZ
  CONTIGUOUS_FSTMT(fluxBufX)
  CONTIGUOUS_FSTMT(fluxBufY)
  CONTIGUOUS_FSTMT(fluxBufZ)
  logical, intent(in), OPTIONAL :: add
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

  real,pointer, dimension(:,:,:,:) :: fluxx,fluxy,fluxz
  CONTIGUOUS_FSTMT(fluxx)
  CONTIGUOUS_FSTMT(fluxy)
  CONTIGUOUS_FSTMT(fluxz)

  logical,parameter :: saveAsFine   = .TRUE.
  logical,parameter :: saveAsCoarse = .TRUE.
  logical,parameter :: saveAsAll    = .TRUE. !save as much locations as we can
  type(Grid_tile_t) :: tileDesc
  integer :: blockID
  integer :: level
  integer :: presVar
  integer :: sx,ex,sy,ey,sz,ez
  logical :: xtrue,ytrue,ztrue
  logical :: multFluxx,multFluxy,multFluxz !whether to multiply by area
  logical :: multFixed !whether to multiply by 2**(NDIM-1)
#if NDIM == 1
  real,parameter :: scalA = 1.0    ! 2.0**(NDIM-1)
#endif
#if NDIM == 2
  real,parameter :: scalA = 2.0    ! 2.0**(NDIM-1)
#endif
#if NDIM == 3
  real,parameter :: scalA = 4.0    ! 2.0**(NDIM-1)
#endif
  integer, dimension(MDIM) :: offs
  integer                  :: offx, offy, offz
  integer                  :: blkLev
  logical                  :: loCase2, hiCase2
  logical                  :: loCase4, hiCase4
  logical                  :: loCase24, hiCase24
  logical                  :: allIsFlux
  real,allocatable,target :: faceAreas(:,:,:)
  real,pointer            :: areaLeft(:,:,:)
  real                    :: wt

  level = blockDesc % level

  wt=0.0
  if (present(add)) then
     if(add)wt=1.0
  end if

  if (present(isFluxDensity)) then
     allIsFlux = ( size(isFluxDensity,1) > 0 .AND. .NOT. ANY(isFluxDensity) )
  else if (present(add)) then
     allIsFlux = .TRUE.
  else
     allIsFlux = .FALSE.
  end if

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

!!$  if(present(axis))then
!!$     xtrue = (axis==IAXIS)
!!$     ytrue = (axis==JAXIS)
!!$     ztrue = (axis==KAXIS)
!!$  end if

  if (allIsFlux) then
     multFluxx = .false. ; multFluxy = .false. ; multFluxz = .false.
     multFixed = .FALSE.
  else
     select case (gr_geometry)
     case(CARTESIAN)
        multFluxx = .false. ; multFluxy = .false. ; multFluxz = .false.
     case(SPHERICAL)
        multFluxx = (NDIM>1); multFluxy = .TRUE.  ; multFluxz = .TRUE.
     case(POLAR)
        multFluxx = .FALSE. ; multFluxy = .FALSE. ; multFluxz = .TRUE.
     case(CYLINDRICAL)
        multFluxx = .FALSE. ; multFluxy = .TRUE.  ; multFluxz = .FALSE.
     end select
     multFixed = ((NDIM > 1) .AND. (gr_geometry .NE. CARTESIAN))
  end if

  if ((level < lrefine_max .OR. saveAsAll .OR. saveAsFine  ) .OR. &
      (level > 1           .OR. saveAsAll .OR. saveAsCoarse)) then
     tileDesc = blockDesc       !We shall, however, assume that the "tile" is a full block. For now.
     blockID=tileDesc%id
     blkLev =tileDesc%level
     fluxx(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxBufX ! fluxx,fluxy,fluxz use local (Paramesh) index counting
     fluxy(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxBufY ! fluxBuf[XYZ] use the global index convention (for the level)
     fluxz(1:,gr_iloFl:,gr_jloFl:,gr_kloFl:) => fluxBufZ
     offs(:) = tileDesc%blkLimitsGC(LOW,1:MDIM) - 1
     offx = offs(IAXIS); offy = offs(JAXIS); offz = offs(KAXIS)

     if(xtrue) then
        flux_x(:nfluxes,1,:,:,blockID) = fluxx(:,sx,sy:ey,sz:ez) + wt*flux_x(:nfluxes,1,:,:,blockID)
        flux_x(:nfluxes,2,:,:,blockID) = fluxx(:,ex+1,sy:ey,sz:ez) + wt*flux_x(:nfluxes,2,:,:,blockID)
        gr_xflx(:,1,:,:,blockID) = fluxx(:,sx+1,sy:ey,sz:ez) + wt*gr_xflx(:,1,:,:,blockID)
        gr_xflx(:,2,:,:,blockID) = fluxx(:,ex,sy:ey,sz:ez)+wt*gr_xflx(:,2,:,:,blockID)
#ifdef FLASH_HYDRO_UNSPLIT
        !! Store transverse components for the faces in global scratch arrays.
#if NDIM >=2
        gr_xflx_yface(:,3:NXB-1, 1     ,:,blockID) = fluxx(:,sx+2:ex-1,sy,sz:ez) +wt*&
             gr_xflx_yface(:,3:NXB-1, 1     ,:,blockID)
        gr_xflx_yface(:,3:NXB-1, 2     ,:,blockID) = fluxx(:,sx+2:ex-1,ey,sz:ez)+ wt*&
             gr_xflx_yface(:,3:NXB-1, 2     ,:,blockID)
#if NDIM>2
        gr_xflx_zface(:,3:NXB-1,2:NYB-1,1,blockID) = fluxx(:,sx+2:ex-1,sy+1:ey-1,sz) +wt*&
             gr_xflx_zface(:,3:NXB-1,2:NYB-1,1,blockID)
        gr_xflx_zface(:,3:NXB-1,2:NYB-1,2,blockID) = fluxx(:,sx+2:ex-1,sy+1:ey-1,ez) +wt*&
             gr_xflx_zface(:,3:NXB-1,2:NYB-1,2,blockID)
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

        ! * We copy to flux_{x,y,z} in cases i,ii,iii,iv (i.e., all the above)
        ! * We multiply with face areas if
        !   - other conditions are satisfied (geometry and direction), and
        !   - face area is nonzero, and
        !   - case is ii or iv.

        if (multFluxx .OR. multFixed) then
           loCase2 = (surr_blks(1,1,1+K2D,1+K3D,blockID) == -1)
           hiCase2 = (surr_blks(1,3,1+K2D,1+K3D,blockID) == -1)
           loCase4 = (surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
                      surr_blks(3,1,1+K2D,1+K3D,blockID) == PARENT_BLK)
           hiCase4 = (surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
                      surr_blks(3,3,1+K2D,1+K3D,blockID) == PARENT_BLK)
           loCase24 = loCase2 .OR. loCase4
           hiCase24 = hiCase2 .OR. hiCase4

           if (multFluxx .AND. (loCase24 .OR. hiCase24)) then
              allocate(faceAreas(offx+sx:offx+ex+1,offy+sy:offy+ey,offz+sz:offz+ez))
              areaLeft(sx:,sy:,sz:)  => faceAreas

              call Grid_getCellFaceAreas(IAXIS, blkLev, &
                                         lbound(faceAreas), ubound(faceAreas), &
                                         faceAreas)

              do presVar = 1,nfluxes
                 if (loCase24) then
                    where (areaLeft(sx,sy:ey,sz:ez).NE.0.0)
                       flux_x(presVar,1,:,:,blockID) = flux_x(presVar,1,:,:,blockID) * areaLeft(sx,sy:ey,sz:ez)
                    end where
                 end if
                 if (hiCase24) then
                    flux_x(presVar,2,:,:,blockID) = flux_x(presVar,2,:,:,blockID) * areaLeft(ex+1,sy:ey,sz:ez)
                 end if
              end do
              deallocate(faceAreas)
           else if (loCase4 .OR. hiCase4) then
              do presVar = 1,nfluxes
                 if (loCase4) then
                    flux_x(presVar,1,:,:,blockID) = flux_x(presVar,1,:,:,blockID) * scalA
                 end if
                 if (hiCase4) then
                    flux_x(presVar,2,:,:,blockID) = flux_x(presVar,2,:,:,blockID) * scalA
                 end if
              end do
           end if
        end if

     end if

#if NDIM>1
     if(ytrue) then
        flux_y(:nfluxes,:,1,:,blockID)  = fluxy(:,sx:ex,sy,sz:ez)  +wt*&
             flux_y(:nfluxes,:,1,:,blockID)
        flux_y(:nfluxes,:,2,:,blockID)  = fluxy(:,sx:ex,ey+1,sz:ez) +wt*&
             flux_y(:nfluxes,:,2,:,blockID)
        gr_yflx(:,:,1,:,blockID) =  fluxy(:,sx:ex,sy+1,sz:ez) +wt*&
             gr_yflx(:,:,1,:,blockID)
        gr_yflx(:,:,2,:,blockID) =  fluxy(:,sx:ex,ey,sz:ez) +wt*&
             gr_yflx(:,:,2,:,blockID)
#ifdef FLASH_HYDRO_UNSPLIT
        !! Store transverse components for the faces in global scratch arrays.
        gr_yflx_xface(:,1,3:NYB-1,:,blockID) = fluxy(:,sx,sy+2:ey-1,sz:ez) +wt*&
             gr_yflx_xface(:,1,3:NYB-1,:,blockID)
        gr_yflx_xface(:,2,3:NYB-1,:,blockID) = fluxy(:,ex,sy+2:ey-1,sz:ez) +wt*&
             gr_yflx_xface(:,2,3:NYB-1,:,blockID)
#if NDIM>2
        gr_yflx_zface(:,2:NXB-1,3:NYB-1,1,blockID) = fluxy(:,sx+1:ex-1,sy+2:ey-1,sz) +wt*&
             gr_yflx_zface(:,2:NXB-1,3:NYB-1,1,blockID)
        gr_yflx_zface(:,2:NXB-1,3:NYB-1,2,blockID) = fluxy(:,sx+1:ex-1,sy+2:ey-1,ez) +wt*&
             gr_yflx_zface(:,2:NXB-1,3:NYB-1,2,blockID)
#endif
#endif
        if (multFluxy .OR. multFixed) then
           loCase2 = (surr_blks(1,2,1,1+K3D,blockID) == -1)
           hiCase2 = (surr_blks(1,2,3,1+K3D,blockID) == -1)
           loCase4 = (surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
                      surr_blks(3,2,1,1+K3D,blockID) == PARENT_BLK)
           hiCase4 = (surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
                      surr_blks(3,2,3,1+K3D,blockID) == PARENT_BLK)
           loCase24 = loCase2 .OR. loCase4
           hiCase24 = hiCase2 .OR. hiCase4

           if (multFluxy .AND. (loCase24 .OR. hiCase24)) then
              allocate(faceAreas(offx+sx:offx+ex,offy+sy:offy+ey+1,offz+sz:offz+ez))
              areaLeft(sx:,sy:,sz:)  => faceAreas

              call Grid_getCellFaceAreas(JAXIS, blkLev, &
                                         lbound(faceAreas), ubound(faceAreas), &
                                         faceAreas)

              do presVar = 1,nfluxes
                 if (loCase24) then
                    where (areaLeft(sx:ex,sy,sz:ez).NE.0.0)
                       flux_y(presVar,:,1,:,blockID) = flux_y(presVar,:,1,:,blockID) * areaLeft(sx:ex,sy,sz:ez)
                    end where
                 end if
                 if (hiCase24) then
                    where (areaLeft(sx:ex,ey+1,sz:ez).NE.0.0)
                       flux_y(presVar,:,2,:,blockID) = flux_y(presVar,:,2,:,blockID) * areaLeft(sx:ex,ey+1,sz:ez)

                    end where
                 end if
              end do
              deallocate(faceAreas)
           else if (loCase4 .OR. hiCase4) then
              do presVar = 1,nfluxes
                 if (loCase4) then
                    flux_y(presVar,:,1,:,blockID) = flux_y(presVar,:,1,:,blockID) * scalA
                 end if
                 if (hiCase4) then
                    flux_y(presVar,:,2,:,blockID) = flux_y(presVar,:,2,:,blockID) * scalA
                 end if
              end do
           end if
        end if
     end if
#endif

#if NDIM>2
     if(ztrue) then
        flux_z(:nfluxes,:,:,1,blockID) = fluxz(:,sx:ex,sy:ey,sz) +wt*&
             flux_z(:nfluxes,:,:,1,blockID)
        flux_z(:nfluxes,:,:,2,blockID) = fluxz(:,sx:ex,sy:ey,ez+1) +wt*&
             flux_z(:nfluxes,:,:,2,blockID)
        gr_zflx(:,:,:,1,blockID) = fluxz(:,sx:ex,sy:ey,sz+1) +wt*&
             gr_zflx(:,:,:,1,blockID)
        gr_zflx(:,:,:,2,blockID) = fluxz(:,sx:ex,sy:ey,ez) +wt*&
             gr_zflx(:,:,:,2,blockID)
#ifdef FLASH_HYDRO_UNSPLIT
        !! Store transverse components for the faces in global scratch arrays.
        gr_zflx_xface(:, 1,      : ,3:NZB-1,blockID) = fluxz(:,sx,sy:ey,sz+2:ez-1) +wt*&
             gr_zflx_xface(:, 1,      : ,3:NZB-1,blockID)
        gr_zflx_xface(:, 2,      : ,3:NZB-1,blockID) = fluxz(:,ex,sy:ey,sz+2:ez-1) +wt*&
             gr_zflx_xface(:, 2,      : ,3:NZB-1,blockID)

        gr_zflx_yface(:,2:NXB-1, 1 ,3:NZB-1,blockID) = fluxz(:,sx+1:ex-1,sy,sz+2:ez-1) +wt*&
             gr_zflx_yface(:,2:NXB-1, 1 ,3:NZB-1,blockID)
        gr_zflx_yface(:,2:NXB-1, 2 ,3:NZB-1,blockID) = fluxz(:,sx+1:ex-1,ey,sz+2:ez-1) +wt*&
             gr_zflx_yface(:,2:NXB-1, 2 ,3:NZB-1,blockID)
#endif
        if (multFluxz .OR. multFixed) then
           loCase2 = (surr_blks(1,2,2,1,blockID) == -1)
           hiCase2 = (surr_blks(1,2,2,3,blockID) == -1)
           loCase4 = (surr_blks(1,2,2,1,blockID) > 0 .AND. &
                      surr_blks(3,2,2,1,blockID) == PARENT_BLK)
           hiCase4 = (surr_blks(1,2,2,3,blockID) > 0 .AND. &
                      surr_blks(3,2,2,3,blockID) == PARENT_BLK)
           loCase24 = loCase2 .OR. loCase4
           hiCase24 = hiCase2 .OR. hiCase4

           if (multFluxz .AND. (loCase24 .OR. hiCase24)) then
              allocate(faceAreas(offx+sx:offx+ex,offy+sy:offy+ey,offz+sz:offz+ez+1))
              areaLeft(sx:,sy:,sz:)  => faceAreas

              call Grid_getCellFaceAreas(KAXIS, blkLev, &
                                         lbound(faceAreas), ubound(faceAreas), &
                                         faceAreas)
              do presVar = 1,nfluxes
                 if (loCase24) then
                    flux_z(presVar,:,:,1,blockID) = flux_z(presVar,:,:,1,blockID) * areaLeft(sx:ex,sy:ey,sz)
                 end if
                 if (hiCase24) then
                    flux_z(presVar,:,:,2,blockID) = flux_z(presVar,:,:,2,blockID) * areaLeft(sx:ex,sy:ey,ez+1)
                 end if
              end do
              deallocate(faceAreas)
           else if (loCase4 .OR. hiCase4) then
              do presVar = 1,nfluxes
                 if (loCase4) then
                    flux_z(presVar,:,:,1,blockID) = flux_z(presVar,:,:,1,blockID) * scalA
                 end if
                 if (hiCase4) then
                    flux_z(presVar,:,:,2,blockID) = flux_z(presVar,:,:,2,blockID) * scalA
                 end if
              end do
           end if
        end if
     end if
#endif
  end if

end subroutine Grid_putFluxData_block
