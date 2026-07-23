!!****if* source/Grid/GridMain/AMR/Amrex/Grid_getFluxCorrData_block
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
!!  Grid_getFluxCorrData_block
!!
!! SYNOPSIS
!!
!!  call Grid_getFluxCorrData_block(type(Grid_tile_t)(in) :: blockDesc,
!!                              real(OUT),TARGET      :: fluxBufX(lo(1): ,lo(2): ,lo(3): ,: ),
!!                              real(OUT),TARGET      :: fluxBufY(lo(1): ,lo(2): ,lo(3): ,: ),
!!                              real(OUT),TARGET      :: fluxBufZ(lo(1): ,lo(2): ,lo(3): ,: ),
!!                              integer(in)           :: lo(3),
!!                              logical(IN), OPTIONAL :: isFluxDensity(:))
!!
!! DESCRIPTION
!!
!!   Get flux corrections from semipermanent flux storage (SPFS).
!!
!!    fluxBuf := "communicated fine fluxes" - "saved coarse fluxes"  AT coarse side of f/c bdry;
!!            :=  0.0                                                ELSEWHERE.
!!
!!   Flux corrections are returned in fluxBuf[XYZ] arguments.
!!
!!   Only fluxes at locations that represent the coarse side of fine/coarse
!!   block boundaries can hold nonzero flux correction data on return.
!!   Other elements of the flux correction buffers fluxBuf[XYZ] are set to zero.
!!
!! ARGUMENTS
!!
!!   blockDesc : describes the current block.
!!               Note that this should be a full block, not a tile representing
!!               a partial block.
!!
!!   fluxBufX :  buffer for flux correction (difference) in IAXIS-direction
!!
!!   fluxBufY :  buffer for flux correction (difference) in JAXIS-direction; output is 0 if NDIM < 2
!!
!!   fluxBufZ :  buffer for flux correction (difference) in KAXIS-direction; output is 0 if NDIM < 3
!!
!!   lo :        lower bounds for the spatial indices of the flux buffers
!!
!!   isFluxDensity : indicates, for each flux component, whether the component
!!                   is a flux proper (if TRUE) or a flux density (otherwise).
!!                   This may be either removed, or changed into a scalar flag,
!!                   later.
!!                   ONLY PARTIALLY SUPPORTED IN THIS IMPLEMENTATION, WILL ABORT
!!                   IF A MIXTURE OF fluxes AND flux density IS INDICATED.
!!
!! NOTES
!!
!!   The arrays fluxBufX, fluxBufY, fluxBufZ are declared with the
!!   index order appropriate for use with AMReX, and are thus not
!!   subject to index reordering.
!!
!!   flux buffer arrays should contain space for fluxes of all valid cells
!!   in the block, excluding guard cells.
!!
!!   This interface does not require level-wide fluxes to be allocated.
!!
!!   SPFS means semi-permanent flux storage. When using a Grid
!!   implementation based on AMReX, SPFS is implemented by an AMReX
!!   flux register class, such as FlashFluxRegister, possibly in
!!   combination with some auxiliary store in module gr_auxFluxData.
!!
!! SEE ALSO
!!
!!   Grid_putFluxData_block
!!   Grid_communicateFluxes
!!   Grid_correctFluxData
!!   Hydro
!!***

#include "Simulation.h"

subroutine Grid_getFluxCorrData_block(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, isFluxDensity)
  use Driver_interface, ONLY : Driver_abort
  use Grid_interface, ONLY : Grid_getCellFaceAreas

  use Grid_tile, ONLY : Grid_tile_t
  use Grid_data,      ONLY : gr_geometry

  use gr_leafBlockInfo, ONLY: gr_getLeafBlockNo
  use gr_auxFluxData, ONLY : gr_iloFl, gr_jloFl, gr_kloFl
  use gr_auxFluxData, ONLY : tflux_x=>gr_tfluxX, tflux_y=>gr_tfluxY, tflux_z=>gr_tfluxZ
  use gr_physicalMultifabs, ONLY : unk, flux_registers
  implicit none

#include "constants.h"

  type(Grid_tile_t), intent(in) :: blockDesc
  integer,intent(in) :: lo(3)
  real,intent(OUT),dimension(lo(1): ,lo(2): ,lo(3): ,: ),TARGET :: fluxBufX, fluxBufY, fluxBufZ
  logical, intent(IN), OPTIONAL :: isFluxDensity(:)

  real,pointer, dimension(:,:,:,:) :: fluxx,fluxy,fluxz
  integer :: coarseLev          !level of this (coarse) block, FLASH convention (1-based)
  integer :: ilev               !level of FlasFluxRegister from whose coarse side we shall load, AMReX convention
  integer :: igrd               !grid index as used by AMReX

  integer :: blockID            !actually, a leaf ID
  integer :: level              ! coarse level, for a block in which fluxes get corrected
  integer :: presVar
  integer :: sx,ex,sy,ey,sz,ez
  logical :: divideFluxx,divideFluxy,divideFluxz !whether to divide by area
  real,allocatable :: faceAreas(:,:,:)
  logical                  :: allIsFlux

#ifndef USE_AMREX_FLASHFLUXREGISTER
  call Driver_abort("Grid_getFluxCorrData_block.F90 requires amrex_flash_fluxregister,&
       & make sure USE_AMREX_FLASHFLUXREGISTER is defined!")
#else

  coarseLev = blockDesc % level
  ilev = coarseLev
  if (ilev >= size(unk,1)) then
     ! There are no corrections available for this level.
     ! The caller could have avoided calling us, but we shall be reasonable.
     fluxBufX(:,:,:,:) = 0.0
     fluxBufY(:,:,:,:) = 0.0
     fluxBufZ(:,:,:,:) = 0.0
     ! NOTHING more to do here!
     ! RETURN imediately!
     return
  end if

  if (present(isFluxDensity)) then
     allIsFlux = ( size(isFluxDensity,1) > 0 .AND. .NOT. ANY(isFluxDensity) )
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

  if (allIsFlux) then
     divideFluxx = .false. ; divideFluxy = .false. ; divideFluxz = .false.
  else
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
  end if

  if (ilev > 0) then

     if (present(isFluxDensity)) then
        if (.NOT.(ALL(isFluxDensity) .OR. allIsFlux)) &
             call Driver_abort("Grid_putFluxData_block: isFluxDensity not yet fully supported")
     end if

     igrd = blockDesc % grid_index

     fluxBufX = 0.0
     if (divideFluxx) then
        allocate(faceAreas(lbound(fluxBufX,1):ubound(fluxBufX,1), &
                           lbound(fluxBufX,2):ubound(fluxBufX,2), &
                           lbound(fluxBufX,3):ubound(fluxBufX,3)))
        call Grid_getCellFaceAreas     (IAXIS,   coarseLev, lbound(fluxBufX),  ubound(fluxBufX),  faceAreas)
        call flux_registers(ilev)%load(fluxBufX,faceAreas, lbound(fluxBufX)-1, ubound(fluxBufX)-1, igrd, 0)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%load(fluxBufX,           lbound(fluxBufX)-1, ubound(fluxBufX)-1, igrd, 0)
     end if
     blockID = gr_getLeafBlockNo(blockDesc)
     fluxx(gr_iloFl:,gr_jloFl:,gr_kloFl:,1:) => fluxBufX ! fluxBuf[XYZ] use the global index convention (for the level)
     if (ANY(fluxx(sx,sy:ey,sz:ez,:) .NE. 0.0)) then
        fluxx(sx,sy:ey,sz:ez,:) = fluxx(sx,sy:ey,sz:ez,:) - tflux_x(1,:,:,:NFLUXES,blockID)
     end if
     if (ANY(fluxx(ex+1,sy:ey,sz:ez,:) .NE. 0.0)) &
          fluxx(ex+1,sy:ey,sz:ez,:) = fluxx(ex+1,sy:ey,sz:ez,:) - tflux_x(2,:,:,:NFLUXES,blockID)

     fluxBufY = 0.0
#if NDIM > 1
     if (divideFluxy) then
        allocate(faceAreas(lbound(fluxBufY,1):ubound(fluxBufY,1), &
                           lbound(fluxBufY,2):ubound(fluxBufY,2), &
                           lbound(fluxBufY,3):ubound(fluxBufY,3)))
        call Grid_getCellFaceAreas     (JAXIS,   coarseLev, lbound(fluxBufY),  ubound(fluxBufY),  faceAreas)
        call flux_registers(ilev)%load(fluxBufY,faceAreas, lbound(fluxBufY)-1, ubound(fluxBufY)-1, igrd, 1)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%load(fluxBufY,           lbound(fluxBufY)-1, ubound(fluxBufY)-1, igrd, 1)
     end if
     fluxy(gr_iloFl:,gr_jloFl:,gr_kloFl:,1:) => fluxBufY ! fluxBuf[XYZ] use the global index convention (for the level)
     if (ANY(fluxy(sx:ex,sy,sz:ez,:) .NE. 0.0)) &
          fluxy(sx:ex,sy,sz:ez,:) = fluxy(sx:ex,sy,sz:ez,:) - tflux_y(:,1,:,:NFLUXES,blockID)
     if (ANY(fluxy(sx:ex,ey+1,sz:ez,:) .NE. 0.0)) &
          fluxy(sx:ex,ey+1,sz:ez,:) = fluxy(sx:ex,ey+1,sz:ez,:) - tflux_y(:,2,:,:NFLUXES,blockID)
#endif

     fluxBufZ = 0.0
#if NDIM == 3
     if (divideFluxz) then
        allocate(faceAreas(lbound(fluxBufZ,1):ubound(fluxBufZ,1), &
                           lbound(fluxBufZ,2):ubound(fluxBufZ,2), &
                           lbound(fluxBufZ,3):ubound(fluxBufZ,3)))
        call Grid_getCellFaceAreas     (KAXIS,   coarseLev, lbound(fluxBufZ),  ubound(fluxBufZ),  faceAreas)
        call flux_registers(ilev)%load(fluxBufZ,faceAreas, lbound(fluxBufZ)-1, ubound(fluxBufZ)-1, igrd, 2)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%load(fluxBufZ,           lbound(fluxBufZ)-1, ubound(fluxBufZ)-1, igrd, 2)
     end if
     fluxz(gr_iloFl:,gr_jloFl:,gr_kloFl:,1:) => fluxBufZ ! fluxBuf[XYZ] use the global index convention (for the level)
     if (ANY(fluxz(sx:ex,sy:ey,sz,:) .NE. 0.0)) &
          fluxz(sx:ex,sy:ey,sz,:) = fluxz(sx:ex,sy:ey,sz,:) - tflux_z(:,:,1,:NFLUXES,blockID)
     if (ANY(fluxz(sx:ex,sy:ey,ez+1,:) .NE. 0.0)) &
          fluxz(sx:ex,sy:ey,ez+1,:) = fluxz(sx:ex,sy:ey,ez+1,:) - tflux_z(:,:,2,:NFLUXES,blockID)
#endif
  end if

#endif
end subroutine Grid_getFluxCorrData_block
