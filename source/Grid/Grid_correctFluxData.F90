!!****f* source/Grid/Grid_correctFluxData
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
!! SEE ALSO
!!
!!   Grid_putFluxData_block
!!   Grid_communicateFluxes
!!   Grid_correctFluxData_xtra
!!   Hydro
!!
!!***

!!REORDER(4): fluxBuf[XYZ]

subroutine Grid_correctFluxData(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, isFluxDensity)
  use Grid_tile, ONLY : Grid_tile_t
  implicit none

#include "FortranLangFeatures.fh"
  type(Grid_tile_t), intent(in) :: blockDesc
  integer,intent(in) :: lo(3)
  real,intent(INOUT),dimension(: , lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxBufX, fluxBufY, fluxBufZ
  CONTIGUOUS_FSTMT(fluxBufX)
  CONTIGUOUS_FSTMT(fluxBufY)
  CONTIGUOUS_FSTMT(fluxBufZ)
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

end subroutine Grid_correctFluxData
