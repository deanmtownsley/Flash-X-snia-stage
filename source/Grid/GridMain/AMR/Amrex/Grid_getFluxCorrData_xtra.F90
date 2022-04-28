!!****if* source/Grid/GridMain/AMR/Amrex/Grid_getFluxCorrData_xtra
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
!!                                 real(IN), TARGET,CONTIGUOUS :: fluxBufX(lo(1): ,lo(2): ,lo(3): ,: ),
!!                                 real(IN), TARGET,CONTIGUOUS :: fluxBufY(lo(1): ,lo(2): ,lo(3): ,: ),
!!                                 real(IN), TARGET,CONTIGUOUS :: fluxBufZ(lo(1): ,lo(2): ,lo(3): ,: ),
!!                                 integer(in)           :: lo(3),
!!                                 real(OUT),TARGET,CONTIGUOUS :: fluxCorrX(lo(1): ,lo(2): ,lo(3): ,: ),
!!                                 real(OUT),TARGET,CONTIGUOUS :: fluxCorrY(lo(1): ,lo(2): ,lo(3): ,: ),
!!                                 real(OUT),TARGET,CONTIGUOUS :: fluxCorrZ(lo(1): ,lo(2): ,lo(3): ,: ),
!!                                 logical(IN), OPTIONAL :: isFluxDensity)
!!
!! DESCRIPTION
!!
!!   Get flux corrections from semipermanent flux storage (SPFS).
!!
!!   Flux corrections are returned in fluxCorr[XYZ] arguments.
!!   The arguments fluxBuf[XYZ] are input only and represent coarse
!!   fluxes, i.e., the ones for which corrections need to be computed.
!!
!!   Only fluxes at locations that represent the coarse side of fine/coarse
!!   block boundaries can hold nonzero flux correction data on return.
!!   Other elements of the flux correction buffers are set to zero.
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
!!   fluxCorrY : flux correction (difference) for JAXIS direction
!!
!!   fluxCorrZ : flux correction (difference) for KAXIS direction
!!
!!   isFluxDensity : indicates, for each flux component, whether the component
!!                   is a flux proper (if TRUE) or a flux density (otherwise).
!!                   This may be either removed, or changed into a scalar flag,
!!                   later.
!!
!! NOTES
!!
!!   The arrays fluxBufX, fluxBufY, fluxBufZ are declared with the
!!   index order appropriate for use with AMReX, and are thus not
!!   subject to index reordering. The same applies to the arrays
!!   fluxCorrX, fluxCorrY, fluxCorrZ.
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
!! SEE ALSO
!!
!!   Grid_putFluxData_block
!!   Grid_communicateFluxes
!!   Grid_correctFluxData
!!   Hydro
!!***

#include "Simulation.h"

subroutine Grid_getFluxCorrData_xtra(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, fluxCorrX,fluxCorrY,fluxCorrZ, isFluxDensity)
  use Driver_interface, ONLY : Driver_abort
  use Grid_interface, ONLY : Grid_getCellFaceAreas

  use Grid_tile, ONLY : Grid_tile_t
  use Grid_data,      ONLY : gr_geometry
  use gr_physicalMultifabs, ONLY : flux_registers
  implicit none

#include "constants.h"
#include "FortranLangFeatures.fh"

  type(Grid_tile_t), intent(in) :: blockDesc
  integer,intent(in) :: lo(3)
  real,intent(in)   ,dimension(lo(1): ,lo(2): ,lo(3): ,: ),TARGET :: fluxBufX, fluxBufY, fluxBufZ
  real,intent(OUT)  ,dimension(lo(1): ,lo(2): ,lo(3): ,: ),TARGET :: fluxCorrX, fluxCorrY, fluxCorrZ
  CONTIGUOUS_FSTMT(fluxBufX)
  CONTIGUOUS_FSTMT(fluxBufY)
  CONTIGUOUS_FSTMT(fluxBufZ)
  CONTIGUOUS_FSTMT(fluxCorrX)
  CONTIGUOUS_FSTMT(fluxCorrY)
  CONTIGUOUS_FSTMT(fluxCorrZ)
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

  real,parameter    :: scaleF = 1.0, scaleC = -1.0
  integer :: coarseLev          !level of this (coarse) block, FLASH convention (1-based)
  integer :: ilev               !level of FlasFluxRegister from whose coarse side we shall load, AMReX convention
  integer :: igrd               !grid index as used by AMReX

  logical :: divideFluxx,divideFluxy,divideFluxz !whether to divide by area
  real,allocatable :: faceAreas(:,:,:)

#ifndef USE_AMREX_FLASHFLUXREGISTER
  call Driver_abort("Grid_getFluxCorrData_xtra.F90 requires amrex_flash_fluxregister,&
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
  if (ilev > 0) then

     if (present(isFluxDensity)) then
        if (.NOT.ALL(isFluxDensity)) &
             call Driver_abort("Grid_getFluxCorrData_xtra: isFluxDensity is not yet supported")
     end if

     igrd = blockDesc % grid_index

     fluxCorrX = 0.0
     if (divideFluxx) then
        allocate(faceAreas(lbound(fluxBufX,1):ubound(fluxBufX,1), &
                           lbound(fluxBufX,2):ubound(fluxBufX,2), &
                           lbound(fluxBufX,3):ubound(fluxBufX,3)))
        call Grid_getCellFaceAreas     (IAXIS,   coarseLev, lbound(fluxBufX),  ubound(fluxBufX),  faceAreas)
        call flux_registers(ilev)%load(fluxCorrX,faceAreas, lbound(fluxBufX)-1,ubound(fluxBufX)-1,fluxBufX, igrd,0, scaleF,scaleC)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%load(fluxCorrX,           lbound(fluxBufX)-1,ubound(fluxBufX)-1,fluxBufX, igrd,0, scaleF,scaleC)
     end if
     fluxCorrY = 0.0
#if NDIM > 1
     if (divideFluxy) then
        allocate(faceAreas(lbound(fluxBufY,1):ubound(fluxBufY,1), &
                           lbound(fluxBufY,2):ubound(fluxBufY,2), &
                           lbound(fluxBufY,3):ubound(fluxBufY,3)))
        call Grid_getCellFaceAreas     (JAXIS,   coarseLev, lbound(fluxBufY),  ubound(fluxBufY),  faceAreas)
        call flux_registers(ilev)%load(fluxCorrY,faceAreas, lbound(fluxBufY)-1,ubound(fluxBufY)-1,fluxBufY, igrd,1, scaleF,scaleC)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%load(fluxCorrY,           lbound(fluxBufY)-1,ubound(fluxBufY)-1,fluxBufY, igrd,1, scaleF,scaleC)
     end if
#endif
     fluxCorrZ = 0.0
#if NDIM == 3
     if (divideFluxz) then
        allocate(faceAreas(lbound(fluxBufZ,1):ubound(fluxBufZ,1), &
                           lbound(fluxBufZ,2):ubound(fluxBufZ,2), &
                           lbound(fluxBufZ,3):ubound(fluxBufZ,3)))
        call Grid_getCellFaceAreas     (KAXIS,   coarseLev, lbound(fluxBufZ),  ubound(fluxBufZ),  faceAreas)
        call flux_registers(ilev)%load(fluxCorrZ,faceAreas, lbound(fluxBufZ)-1,ubound(fluxBufZ)-1,fluxBufZ, igrd,2, scaleF,scaleC)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%load(fluxCorrZ,           lbound(fluxBufZ)-1,ubound(fluxBufZ)-1,fluxBufZ, igrd,2, scaleF,scaleC)
     end if
#endif
  end if

#endif
end subroutine Grid_getFluxCorrData_xtra
