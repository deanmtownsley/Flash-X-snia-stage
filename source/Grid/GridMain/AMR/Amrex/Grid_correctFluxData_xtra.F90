!!****if* source/Grid/GridMain/AMR/Amrex/Grid_correctFluxData_xtra
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
!!                            real(in),   dimension(lo(1): ,lo(2): ,lo(3): ,:),TARGET,CONTIGUOUS :: fluxOldX,
!!                            real(in),   dimension(lo(1): ,lo(2): ,lo(3): ,:),TARGET,CONTIGUOUS :: fluxOldY,
!!                            real(in),   dimension(lo(1): ,lo(2): ,lo(3): ,:),TARGET,CONTIGUOUS :: fluxOldZ,
!!                            logical(IN), OPTIONAL  :: isFluxDensity)
!!
!!  call Grid_correctFluxData(type(Grid_tile_t)(in) :: blockDesc,
!!                            real(IN)              :: scaleF,
!!                            real(INOUT),TARGET :: fluxBufX(size(fluxOldX,1),size(fluxOldX,2),size(fluxOldX,3),size(fluxOldX,4)),
!!                            real(INOUT),TARGET :: fluxBufY(size(fluxOldY,1),size(fluxOldY,2),size(fluxOldY,3),size(fluxOldY,4)),
!!                            real(INOUT),TARGET :: fluxBufZ(size(fluxOldZ,1),size(fluxOldZ,2),size(fluxOldZ,3),size(fluxOldZ,4)),
!!                            integer(in)           :: lo(MDIM),
!!                            real(IN)              :: scaleC,
!!                            real(in),   dimension(lo(1): ,lo(2): ,lo(3): ,:),TARGET,CONTIGUOUS :: fluxOldX,
!!                            real(in),   dimension(lo(1): ,lo(2): ,lo(3): ,:),TARGET,CONTIGUOUS :: fluxOldY,
!!                            real(in),   dimension(lo(1): ,lo(2): ,lo(3): ,:),TARGET,CONTIGUOUS :: fluxOldZ,
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
!!   The arrays fluxBufX, fluxBufY, fluxBufZ are declared with the
!!   index order appropriate for use with AMReX, and are thus not
!!   subject to index reordering. The same applies to the arrays
!!   fluxOldX, fluxOldY, fluxOldZ.
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
!! SEE ALSO
!!
!!   Grid_putFluxData_block
!!   Grid_communicateFluxes
!!   Grid_correctFluxData
!!   Hydro
!!***

subroutine Grid_correctFluxData_xtra(blockDesc, scaleF, fluxBufX,fluxBufY,fluxBufZ, lo, &
                                                scaleC, fluxOldX,fluxOldY,fluxOldZ,     &
                                                isFluxDensity)
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
  real,intent(in)    :: scaleF,scaleC
  real,intent(in   ),dimension(lo(1): ,lo(2): ,lo(3):, :),TARGET :: fluxOldX, fluxOldY, fluxOldZ
  CONTIGUOUS_FSTMT(fluxOldX)
  CONTIGUOUS_FSTMT(fluxOldY)
  CONTIGUOUS_FSTMT(fluxOldZ)
  real,INTENT(INOUT),dimension(size(fluxOldX,1),size(fluxOldX,2),size(fluxOldX,3),size(fluxOldX,4)),TARGET :: fluxBufX
  real,INTENT(INOUT),dimension(size(fluxOldY,1),size(fluxOldY,2),size(fluxOldY,3),size(fluxOldY,4)),TARGET :: fluxBufY
  real,INTENT(INOUT),dimension(size(fluxOldZ,1),size(fluxOldZ,2),size(fluxOldZ,3),size(fluxOldZ,4)),TARGET :: fluxBufZ
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

  integer :: coarseLev          !level of this (coarse) block, FLASH convention (1-based)
  integer :: ilev               !level of FlashFluxRegister from whose coarse side we shall load, AMReX convention
  integer :: igrd               !grid index as used by AMReX

  logical :: divideFluxx,divideFluxy,divideFluxz !whether to divide by area
  real,allocatable :: faceAreas(:,:,:)

#ifndef USE_AMREX_FLASHFLUXREGISTER
  call Driver_abort("Grid_correctFluxData_xtra.F90 requires amrex_flash_fluxregister,&
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
             call Driver_abort("Grid_correctFluxData_xtra: isFluxDensity is not yet supported")
     end if

     igrd = blockDesc % grid_index

     if (divideFluxx) then
        allocate(faceAreas(lbound(fluxOldX,1):ubound(fluxOldX,1), &
                           lbound(fluxOldX,2):ubound(fluxOldX,2), &
                           lbound(fluxOldX,3):ubound(fluxOldX,3)))
        call Grid_getCellFaceAreas     (IAXIS,   coarseLev, lbound(fluxOldX),  ubound(fluxOldX),  faceAreas)
        call flux_registers(ilev)%load(fluxBufX, faceAreas, lbound(fluxOldX)-1,ubound(fluxOldX)-1,fluxOldX, igrd,0, scaleF,scaleC)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%load(fluxBufX,            lbound(fluxOldX)-1,ubound(fluxOldX)-1,fluxOldX, igrd,0, scaleF,scaleC)
     endif
#if NDIM > 1
     if (divideFluxy) then
        allocate(faceAreas(lbound(fluxOldY,1):ubound(fluxOldY,1), &
                           lbound(fluxOldY,2):ubound(fluxOldY,2), &
                           lbound(fluxOldY,3):ubound(fluxOldY,3)))
        call Grid_getCellFaceAreas     (JAXIS,   coarseLev, lbound(fluxOldY),  ubound(fluxOldY),  faceAreas)
        call flux_registers(ilev)%load(fluxBufY, faceAreas, lbound(fluxOldY)-1,ubound(fluxOldY)-1,fluxOldY, igrd,1, scaleF,scaleC)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%load(fluxBufY,            lbound(fluxOldY)-1,ubound(fluxOldY)-1,fluxOldY, igrd,1, scaleF,scaleC)
     endif
#endif
#if NDIM == 3
     if (divideFluxz) then
        allocate(faceAreas(lbound(fluxOldZ,1):ubound(fluxOldZ,1), &
                           lbound(fluxOldZ,2):ubound(fluxOldZ,2), &
                           lbound(fluxOldZ,3):ubound(fluxOldZ,3)))
        call Grid_getCellFaceAreas     (KAXIS,   coarseLev, lbound(fluxOldZ),  ubound(fluxOldZ),  faceAreas)
        call flux_registers(ilev)%load(fluxBufZ, faceAreas, lbound(fluxOldZ)-1,ubound(fluxOldZ)-1,fluxOldZ, igrd,2, scaleF,scaleC)
        deallocate(faceAreas)
     else
        call flux_registers(ilev)%load(fluxBufZ,            lbound(fluxOldZ)-1,ubound(fluxOldZ)-1,fluxOldZ, igrd,2, scaleF,scaleC)
     endif
#endif
  end if

#endif
end subroutine Grid_correctFluxData_xtra
