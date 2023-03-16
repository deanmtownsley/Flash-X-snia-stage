!!****f* source/Grid/Grid_putFluxData_block
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
!!   A stub that does nothing is used for the Uniform Grid implementation.
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
!!   implementation based on AMReX, SPFS is implemented by an AMReX
!!   flux register class, such as FlashFluxRegister. When using a Grid
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
!!***

!!REORDER(4): fluxBuf[XYZ]

subroutine Grid_putFluxData_block(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, add, isFluxDensity)
  use Grid_tile, ONLY : Grid_tile_t
  implicit none
#include "FortranLangFeatures.fh"
  type(Grid_tile_t), intent(in) :: blockDesc
  integer,intent(in) :: lo(3)
  real,intent(in),dimension(:, lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxBufX,fluxBufY,fluxBufZ
  CONTIGUOUS_FSTMT(fluxBufX)
  CONTIGUOUS_FSTMT(fluxBufY)
  CONTIGUOUS_FSTMT(fluxBufZ)
  logical, intent(in), OPTIONAL :: add
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

  return

end subroutine Grid_putFluxData_block
