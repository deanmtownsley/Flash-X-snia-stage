!!****if* source/Grid/GridMain/AMR/Amrex/Grid_correctFluxData
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
!!                            real(INOUT),TARGET,CONTIGUOUS, dimension(lo(1): ,lo(2): ,lo(3):, :) :: fluxBufX,
!!                            real(INOUT),TARGET,CONTIGUOUS, dimension(lo(1): ,lo(2): ,lo(3):, :) :: fluxBufY,
!!                            real(INOUT),TARGET,CONTIGUOUS, dimension(lo(1): ,lo(2): ,lo(3):, :) :: fluxBufZ,
!!                            integer(in) :: lo(MDIM),
!!                            logical(IN), OPTIONAL  :: isFluxDensity)
!!
!! DESCRIPTION
!!
!!   Correct data in flux arrays by replacing fluxes in certain locations with
!!   data from a higher refinement level.
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
!!   The arrays fluxBufX, fluxBufY, fluxBufZ are declared with the
!!   index order appropriate for use with AMReX, and are thus not
!!   subject to index reordering.
!!
!!   This interface does not require level-wide fluxes to be allocated.
!!
!!   SPFS means semi-permanent flux storage. When using a Grid
!!   implementation based on AMReX, SPFS is implemented by an AMReX
!!   flux register class, such as FlashFluxRegister.
!!
!!   This implementation uses the FlashFluxRegisters implementation of AMReX,
!!   introduced there in Jan/Feb 2020, as underlying flux register implementation.
!!
!! SEE ALSO
!!
!!   Grid_putFluxData_block
!!   Grid_communicateFluxes
!!   Hydro
!!
!!***

subroutine Grid_correctFluxData(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, isFluxDensity)
  use Driver_interface, ONLY : Driver_abort
  use Grid_interface, ONLY : Grid_getCellFaceAreas

  use Grid_tile, ONLY : Grid_tile_t
  use Grid_data,      ONLY : gr_geometry
  use gr_physicalMultifabs, ONLY : flux_registers
  use amrex_amrcore_module, ONLY : amrex_get_finest_level
  implicit none
#include "Simulation.h"
#include "FortranLangFeatures.fh"
#include "constants.h"
  type(Grid_tile_t), intent(in) :: blockDesc
  integer,intent(in) :: lo(3)
  real,intent(INOUT),dimension(lo(1): ,lo(2): ,lo(3): , :),TARGET :: fluxBufX,fluxBufY,fluxBufZ
  CONTIGUOUS_FSTMT(fluxBufX)
  CONTIGUOUS_FSTMT(fluxBufY)
  CONTIGUOUS_FSTMT(fluxBufZ)
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

  integer :: coarseLev          !level of this (coarse) block, FLASH convention (1-based)
  integer :: ilev               !level of FlashFluxRegister from whose coarse side we shall load, AMReX convention
  integer :: igrd               !grid index as used by AMReX

  logical :: divideFluxx,divideFluxy,divideFluxz !whether to divide by area
  real,allocatable :: faceAreas(:,:,:)

#ifndef USE_AMREX_FLASHFLUXREGISTER
  call Driver_abort("Grid_correctFluxData.F90 requires amrex_flash_fluxregister,&
       & make sure USE_AMREX_FLASHFLUXREGISTER is defined!")
#else

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

  coarseLev = blockDesc % level
  ilev = coarseLev
  if (ilev > amrex_get_finest_level())     RETURN

  if (ilev > 0) then

     if (present(isFluxDensity)) then
        if (.NOT.ALL(isFluxDensity)) &
             call Driver_abort("Grid_correctFluxData: isFluxDensity is not yet supported")
     end if

     igrd = blockDesc % grid_index

     if (divideFluxx) then
        allocate(faceAreas(lbound(fluxBufX,1):ubound(fluxBufX,1), &
                           lbound(fluxBufX,2):ubound(fluxBufX,2), &
                           lbound(fluxBufX,3):ubound(fluxBufX,3)))
        call Grid_getCellFaceAreas     (IAXIS,   coarseLev, lbound(fluxBufX),  ubound(fluxBufX),  faceAreas)
        call flux_registers(ilev)%load(fluxBufX, faceAreas, lbound(fluxBufX)-1, ubound(fluxBufX)-1, igrd, 0)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%load(fluxBufX,            lbound(fluxBufX)-1, ubound(fluxBufX)-1, igrd, 0)
     end if
#if NDIM > 1
     if (divideFluxy) then
        allocate(faceAreas(lbound(fluxBufY,1):ubound(fluxBufY,1), &
                           lbound(fluxBufY,2):ubound(fluxBufY,2), &
                           lbound(fluxBufY,3):ubound(fluxBufY,3)))
        call Grid_getCellFaceAreas     (JAXIS,   coarseLev, lbound(fluxBufY),  ubound(fluxBufY),  faceAreas)
        call flux_registers(ilev)%load(fluxBufY, faceAreas, lbound(fluxBufY)-1, ubound(fluxBufY)-1, igrd, 1)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%load(fluxBufY,            lbound(fluxBufY)-1, ubound(fluxBufY)-1, igrd, 1)
     end if
#endif
#if NDIM == 3
     if (divideFluxz) then
        allocate(faceAreas(lbound(fluxBufZ,1):ubound(fluxBufZ,1), &
                           lbound(fluxBufZ,2):ubound(fluxBufZ,2), &
                           lbound(fluxBufZ,3):ubound(fluxBufZ,3)))
        call Grid_getCellFaceAreas     (KAXIS,   coarseLev, lbound(fluxBufZ),  ubound(fluxBufZ),  faceAreas)
        call flux_registers(ilev)%load(fluxBufZ, faceAreas, lbound(fluxBufZ)-1, ubound(fluxBufZ)-1, igrd, 2)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%load(fluxBufZ,            lbound(fluxBufZ)-1, ubound(fluxBufZ)-1, igrd, 2)
     end if
#endif
  end if

#endif
end subroutine Grid_correctFluxData
