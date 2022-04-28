!!****f* source/Grid/Grid_getFluxCorrData_block
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
!!  Grid_getFluxCorrData_block
!!
!! SYNOPSIS
!!
!!  call Grid_getFluxCorrData_block(type(Grid_tile_t)(in) :: blockDesc,
!!                              real(OUT),TARGET      :: fluxBufX(:,lo(1): ,lo(2): ,lo(3): ),
!!                              real(OUT),TARGET      :: fluxBufY(:,lo(1): ,lo(2): ,lo(3): ),
!!                              real(OUT),TARGET      :: fluxBufZ(:,lo(1): ,lo(2): ,lo(3): ),
!!                              integer(in)           :: lo(3),
!!                              logical(IN), OPTIONAL :: isFluxDensity)
!!
!! DESCRIPTION
!!
!!   Get  flux corrections from semipermanent flux storage (SPFS).
!!
!!    fluxBuf := "communicated fine fluxes" - "saved coarse fluxes"  AT coarse side of f/c bdry;
!!            :=  0.0                                                ELSEWHERE.
!!
!!   Only fluxes at locations that represent the coarse side of fine/coarse
!!   block boundaries can hold nonzero flux correction data on return.
!!   Other elements of the flux buffers are set to zero if they represent
!!   faces of any cells that touch a block boundary; data faces of cells
!!   that are farther away from a block boundary are left undefined.
!!
!! ARGUMENTS
!!
!!   blockDesc : describes the current block.
!!               Note that this should be a full block, not a tile representing
!!               a partial block.
!!
!!   fluxBufX :  buffer for fluxes in IAXIS-direction
!!
!!   fluxBufY :  buffer for fluxes in JAXIS-direction;
!!               output should be considered undefined if NDIM < 2.
!!
!!   fluxBufZ :  buffer for fluxes in KAXIS-direction;
!!               output should be considered undefined if NDIM < 3.
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
!!   This interface is currently only implemented for Paramesh4 !
!!
!! SEE ALSO
!!
!!   Grid_putFluxData_block
!!   Grid_communicateFluxes
!!   Grid_correctFluxData
!!   Hydro
!!***

!!REORDER(4): fluxBuf[XYZ]

subroutine Grid_getFluxCorrData_block(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, isFluxDensity)
  use Grid_tile, ONLY : Grid_tile_t

  implicit none

  type(Grid_tile_t), intent(in) :: blockDesc
  integer,intent(in) :: lo(3)
  real,intent(OUT),dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxBufX, fluxBufY, fluxBufZ
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

  fluxBufX = 0.0
  fluxBufY = 0.0
  fluxBufZ = 0.0

end subroutine Grid_getFluxCorrData_block
