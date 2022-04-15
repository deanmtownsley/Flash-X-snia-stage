!!****if* source/Grid/GridMain/AMR/Amrex/Grid_putFluxData_block
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
!!  Grid_putFluxData_block
!!
!! SYNOPSIS
!!
!!  call Grid_putFluxData_block(type(Grid_tile_t)(in) :: blockDesc,
!!                              real(in),Contiguous,TARGET :: fluxBufX(lo(1): ,lo(2): ,lo(3): ,: ),
!!                              real(in),Contiguous,TARGET :: fluxBufY(lo(1): ,lo(2): ,lo(3): ,: ),
!!                              real(in),Contiguous,TARGET :: fluxBufZ(lo(1): ,lo(2): ,lo(3): ,: ),
!!                              integer(in) :: lo(3),
!!                              logical(IN), OPTIONAL :: isFluxDensity)
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
!!   isFluxDensity : indicates, for each flux component, whether the component
!!                   is a flux proper (if TRUE) or a flux density (otherwise).
!!                   This may be either removed, or changed into a scalar flag,
!!                   later.
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

subroutine Grid_putFluxData_block(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo,isFluxDensity)
  use Driver_interface, ONLY : Driver_abort
  use Grid_interface, ONLY : Grid_getCellFaceAreas
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_data,      ONLY : gr_geometry
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
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

  integer :: fineLev            !level of this block, FLASH convention (1-based)
  integer :: ilev               !level of FlasFluxRegister into whose fine side we shall store, AMReX convention
  integer :: igrd               !grid index as used by AMReX

  logical :: multFluxx,multFluxy,multFluxz !whether to multiply by area
  real,allocatable :: faceAreas(:,:,:)

#ifndef USE_AMREX_FLASHFLUXREGISTER
  call Driver_abort("Grid_putFluxData_block.F90 requires amrex_flash_fluxregister,&
       & make sure USE_AMREX_FLASHFLUXREGISTER is defined!")
#else

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

  fineLev = blockDesc % level
  ilev = fineLev - 1
  if (ilev > 0) then

     if (present(isFluxDensity)) then
        if (.NOT.ALL(isFluxDensity)) &
             call Driver_abort("Grid_putFluxData_block: isFluxDensity is not yet supported")
     end if

     igrd = blockDesc % grid_index

     if (multFluxx) then
        allocate(faceAreas(lbound(fluxBufX,1):ubound(fluxBufX,1), &
                           lbound(fluxBufX,2):ubound(fluxBufX,2), &
                           lbound(fluxBufX,3):ubound(fluxBufX,3)))
        call Grid_getCellFaceAreas     (IAXIS,    fineLev,   lbound(fluxBufX),   ubound(fluxBufX), faceAreas)
        call flux_registers(ilev)%store(fluxBufX, faceAreas, lbound(fluxBufX)-1, ubound(fluxBufX)-1, igrd, 0)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%store(fluxBufX,            lbound(fluxBufX)-1, ubound(fluxBufX)-1, igrd, 0)
     end if
#if NDIM > 1
     if (multFluxy) then
        allocate(faceAreas(lbound(fluxBufY,1):ubound(fluxBufY,1), &
                           lbound(fluxBufY,2):ubound(fluxBufY,2), &
                           lbound(fluxBufY,3):ubound(fluxBufY,3)))
        call Grid_getCellFaceAreas     (JAXIS,    fineLev,   lbound(fluxBufY),   ubound(fluxBufY), faceAreas)

        call flux_registers(ilev)%store(fluxBufY, faceAreas, lbound(fluxBufY)-1, ubound(fluxBufY)-1, igrd, 1)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%store(fluxBufY,            lbound(fluxBufY)-1, ubound(fluxBufY)-1, igrd, 1)
     end if
#endif
#if NDIM == 3
     if (multFluxz) then
        allocate(faceAreas(lbound(fluxBufZ,1):ubound(fluxBufZ,1), &
                           lbound(fluxBufZ,2):ubound(fluxBufZ,2), &
                           lbound(fluxBufZ,3):ubound(fluxBufZ,3)))
        call Grid_getCellFaceAreas     (KAXIS,    fineLev,   lbound(fluxBufZ),   ubound(fluxBufZ), faceAreas)
        call flux_registers(ilev)%store(fluxBufZ, faceAreas, lbound(fluxBufZ)-1, ubound(fluxBufZ)-1, igrd, 2)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%store(fluxBufZ,            lbound(fluxBufZ)-1, ubound(fluxBufZ)-1, igrd, 2)
     end if
#endif
  end if

#endif
end subroutine Grid_putFluxData_block
