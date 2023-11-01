!!****if* source/Grid/GridMain/AMR/Amrex/Grid_putFluxData_block
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
!!                              real(in),Contiguous,TARGET :: fluxBufX(lo(1): ,lo(2): ,lo(3): ,: ),
!!                              real(in),Contiguous,TARGET :: fluxBufY(lo(1): ,lo(2): ,lo(3): ,: ),
!!                              real(in),Contiguous,TARGET :: fluxBufZ(lo(1): ,lo(2): ,lo(3): ,: ),
!!                              integer(in) :: lo(3),
!!                              logical(IN), OPTIONAL :: add,
!!                              logical(IN), OPTIONAL :: isFluxDensity(:))
!!
!! DESCRIPTION
!!
!!   Save fluxes passed in as arguments in semipermanent flux storage (SPFS).
!!
!!   For Amrex Grid implementations, SPFS translates to flux registers.
!!
!!   This implementation uses the FlashFluxRegisters implementation of AMReX,
!!   introduced there in Jan/Feb 2020, as underlying flux register implementation.
!!
!! ARGUMENTS
!!
!!   blockDesc : describes the current block.
!!               Note that this should be a full block, not a tile representing
!!               a partial block. For now.
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
!!                   ONLY PARTIALLY SUPPORTED IN THIS IMPLEMENTATION, WILL ABORT
!!                   IF A MIXTURE OF fluxes AND flux density IS INDICATED.
!!
!! NOTES
!!
!!   The arrays fluxBufX, fluxBufY, fluxBufZ are here declared with
!!   the index order appropriate for use with AMReX, and are thus not
!!   subject to index reordering.
!!
!! SEE ALSO
!!
!!   Grid_communicateFluxes
!!   Grid_correctFluxData
!!   Hydro
!!***

subroutine Grid_putFluxData_block(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, add, isFluxDensity)
  use Driver_interface, ONLY : Driver_abort
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Grid_interface, ONLY : Grid_getCellFaceAreas
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_data,      ONLY : gr_geometry
  use gr_leafBlockInfo, ONLY: gr_getLeafBlockNo
  use gr_auxFluxData, ONLY : gr_iloFl, gr_jloFl, gr_kloFl
  use gr_auxFluxData, ONLY : tflux_x=>gr_tfluxX, tflux_y=>gr_tfluxY, tflux_z=>gr_tfluxZ
  use gr_physicalMultifabs, ONLY : flux_registers
  implicit none
#include "Simulation.h"
#include "FortranLangFeatures.fh"
#include "constants.h"
  type(Grid_tile_t), intent(in) :: blockDesc
  integer,intent(in) :: lo(3)
  real,intent(in),dimension(lo(1): ,lo(2): ,lo(3): , :),TARGET :: fluxBufX,fluxBufY,fluxBufZ
  CONTIGUOUS_FSTMT(fluxBufX)
  CONTIGUOUS_FSTMT(fluxBufY)
  CONTIGUOUS_FSTMT(fluxBufZ)
  logical, intent(in), OPTIONAL :: add
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

  real,pointer, dimension(:,:,:,:) :: fluxx,fluxy,fluxz
  integer :: fineLev            !level of this block, FLASH convention (1-based)
  integer :: ilev               !level of FlasFluxRegister into whose fine side we shall store, AMReX convention
  integer :: igrd               !grid index as used by AMReX

  integer :: blockID            !actually, a leaf ID
  integer :: sx,ex,sy,ey,sz,ez
  logical :: multFluxx,multFluxy,multFluxz !whether to multiply by area
#if NDIM == 1
  real,parameter :: scalA = 1.0    ! 2.0**(NDIM-1)
#endif
#if NDIM == 2
  real,parameter :: scalA = 2.0    ! 2.0**(NDIM-1)
#endif
#if NDIM == 3
  real,parameter :: scalA = 4.0    ! 2.0**(NDIM-1)
#endif
  logical                  :: allIsFlux
  real,allocatable :: faceAreas(:,:,:)
  real             :: wt
  logical          :: addit

#ifndef USE_AMREX_FLASHFLUXREGISTER
  call Driver_abort("Grid_putFluxData_block.F90 requires amrex_flash_fluxregister,&
       & make sure USE_AMREX_FLASHFLUXREGISTER is defined!")
#else

  wt=0.0
  addit = .FALSE.
  if (present(add)) then
     if(add)wt=1.0
     addit = add
  end if

#ifdef DEBUG_GRID
  if (addit) then
!!$     call Logfile_stampMessage("Perhaps you should not try to use nontelescoping Spark with AMReX.")
!!$     print*,'Grid_putFluxData_block: LO,ADD,present(isFluxDensity):', lo, add, present(isFluxDensity)
!!$     if (present(isFluxDensity)) print*,'Grid_putFluxData_block: isFluxDensity:', isFluxDensity
!!$     call Driver_abort("Grid_putFluxData_block: adding not properly tested with Amrex Grid!")
  end if
#endif

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

  if (allIsFlux) then
     multFluxx = .false. ; multFluxy = .false. ; multFluxz = .false.
  else
     select case (gr_geometry)
     case(CARTESIAN)
        multFluxx = .false. ; multFluxy = .false. ; multFluxz = .false.
     case(SPHERICAL)
        multFluxx = (NDIM>1); multFluxy = .TRUE.  ; multFluxz = .TRUE.
     case(POLAR)
        multFluxx = .FALSE. ; multFluxy = .FALSE. ; multFluxz = .TRUE.
     case default               !CYLINDRICAL
        multFluxx = .FALSE. ; multFluxy = .TRUE.  ; multFluxz = .FALSE.
     end select
  end if

  fineLev = blockDesc % level
  ilev = fineLev - 1
  if (ilev > 0) then

     if (present(isFluxDensity)) then
        if (.NOT.(ALL(isFluxDensity) .OR. allIsFlux)) &
             call Driver_abort("Grid_putFluxData_block: isFluxDensity not yet fully supported")
     end if

     igrd = blockDesc % grid_index

     blockID = gr_getLeafBlockNo(blockDesc)
     fluxx(gr_iloFl:,gr_jloFl:,gr_kloFl:,1:) => fluxBufX ! fluxBuf[XYZ] use the global index convention (for the level)
     tflux_x(1,:,:,:NFLUXES,blockID) = &
          fluxx(sx,sy:ey,sz:ez,:) + wt*tflux_x(1,:,:,:NFLUXES,blockID)
     tflux_x(2,:,:,:NFLUXES,blockID) = &
          fluxx(ex+1,sy:ey,sz:ez,:) + wt*tflux_x(2,:,:,:NFLUXES,blockID)
     if (multFluxx) then
        allocate(faceAreas(lbound(fluxBufX,1):ubound(fluxBufX,1), &
                           lbound(fluxBufX,2):ubound(fluxBufX,2), &
                           lbound(fluxBufX,3):ubound(fluxBufX,3)))
        call Grid_getCellFaceAreas     (IAXIS,    fineLev,   lbound(fluxBufX),   ubound(fluxBufX), faceAreas)
        call flux_registers(ilev)%store(fluxBufX, faceAreas, lbound(fluxBufX)-1, ubound(fluxBufX)-1, igrd, 0, addit=add)
        deallocate(faceAreas)
     else if (allIsFlux .AND. NDIM > 1) then
        call flux_registers(ilev)%store(fluxBufX,            lbound(fluxBufX)-1, ubound(fluxBufX)-1, igrd, 0, &
             addit=add,scale=scalA)
     else
        call flux_registers(ilev)%store(fluxBufX,            lbound(fluxBufX)-1, ubound(fluxBufX)-1, igrd, 0, addit=add)
     end if
#if NDIM > 1
     fluxy(gr_iloFl:,gr_jloFl:,gr_kloFl:,1:) => fluxBufY ! fluxBuf[XYZ] use the global index convention (for the level)
     tflux_y(:,1,:,:NFLUXES,blockID)  = fluxy(sx:ex,sy,sz:ez,:)  +wt*&
             tflux_y(:,1,:,:NFLUXES,blockID)
     tflux_y(:,2,:,:NFLUXES,blockID)  = fluxy(sx:ex,ey+1,sz:ez,:) +wt*&
             tflux_y(:,2,:,:NFLUXES,blockID)
     if (multFluxy) then
        allocate(faceAreas(lbound(fluxBufY,1):ubound(fluxBufY,1), &
                           lbound(fluxBufY,2):ubound(fluxBufY,2), &
                           lbound(fluxBufY,3):ubound(fluxBufY,3)))
        call Grid_getCellFaceAreas     (JAXIS,    fineLev,   lbound(fluxBufY),   ubound(fluxBufY), faceAreas)

        call flux_registers(ilev)%store(fluxBufY, faceAreas, lbound(fluxBufY)-1, ubound(fluxBufY)-1, igrd, 1, addit=add)
        deallocate(faceAreas)
     else if (allIsFlux) then
        call flux_registers(ilev)%store(fluxBufY,            lbound(fluxBufY)-1, ubound(fluxBufY)-1, igrd, 1, &
             addit=add,scale=scalA)
     else
        call flux_registers(ilev)%store(fluxBufY,            lbound(fluxBufY)-1, ubound(fluxBufY)-1, igrd, 1, addit=add)
     end if
#endif
#if NDIM == 3
     fluxz(gr_iloFl:,gr_jloFl:,gr_kloFl:,1:) => fluxBufZ ! fluxBuf[XYZ] use the global index convention (for the level)
     tflux_z(:,:,1,:NFLUXES,blockID) = fluxz(sx:ex,sy:ey,sz,:) +wt*&
             tflux_z(:,:,1,:NFLUXES,blockID)
     tflux_z(:,:,2,:NFLUXES,blockID) = fluxz(sx:ex,sy:ey,ez+1,:) +wt*&
             tflux_z(:,:,2,:NFLUXES,blockID)
     if (multFluxz) then
        allocate(faceAreas(lbound(fluxBufZ,1):ubound(fluxBufZ,1), &
                           lbound(fluxBufZ,2):ubound(fluxBufZ,2), &
                           lbound(fluxBufZ,3):ubound(fluxBufZ,3)))
        call Grid_getCellFaceAreas     (KAXIS,    fineLev,   lbound(fluxBufZ),   ubound(fluxBufZ), faceAreas)
        call flux_registers(ilev)%store(fluxBufZ, faceAreas, lbound(fluxBufZ)-1, ubound(fluxBufZ)-1, igrd, 2, addit=add)
        deallocate(faceAreas)
     else if (allIsFlux) then
        call flux_registers(ilev)%store(fluxBufZ,            lbound(fluxBufZ)-1, ubound(fluxBufZ)-1, igrd, 2, &
             addit=add,scale=scalA)
     else
        call flux_registers(ilev)%store(fluxBufZ,            lbound(fluxBufZ)-1, ubound(fluxBufZ)-1, igrd, 2, addit=add)
     end if
#endif
  end if

#endif
end subroutine Grid_putFluxData_block
