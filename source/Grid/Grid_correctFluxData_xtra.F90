!!****f* source/Grid/Grid_correctFluxData_xtra
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
!!                            real(in),   TARGET,CONTIGUOUS :: fluxOldX (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            real(in),   TARGET,CONTIGUOUS :: fluxOldY (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            real(in),   TARGET,CONTIGUOUS :: fluxOldZ (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            logical(IN), OPTIONAL  :: isFluxDensity)
!!
!!  call Grid_correctFluxData(type(Grid_tile_t)(in) :: blockDesc,
!!                            real(IN)              :: scaleF,
!!                            real(INOUT),TARGET :: fluxBufX(size(fluxOldX,1),size(fluxOldX,2),size(fluxOldX,3),size(fluxOldX,4)),
!!                            real(INOUT),TARGET :: fluxBufY(size(fluxOldY,1),size(fluxOldY,2),size(fluxOldY,3),size(fluxOldY,4)),
!!                            real(INOUT),TARGET :: fluxBufZ(size(fluxOldZ,1),size(fluxOldZ,2),size(fluxOldZ,3),size(fluxOldZ,4)),
!!                            integer(in)           :: lo(MDIM),
!!                            real(IN)              :: scaleC,
!!                            real(in),   TARGET,CONTIGUOUS :: fluxOldX (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            real(in),   TARGET,CONTIGUOUS :: fluxOldY (: ,lo(1): ,lo(2): ,lo(3): ),
!!                            real(in),   TARGET,CONTIGUOUS :: fluxOldZ (: ,lo(1): ,lo(2): ,lo(3): ),
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
!!   The arrays fluxOldX, fluxOldY, fluxOldZ are subject to index reordering.
!!   The arrays fluxBufX, fluxBufY, fluxBufZ are index-reorder indirectly,
!!   by having bounds that depend  on the shape on the corresponding
!!   fluxOldX, fluxOldY, fluxOldZ arrays.
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
!!
!!
!!***

!!REORDER(4): fluxOld[XYZ]

subroutine Grid_correctFluxData_xtra(blockDesc, scaleF, fluxBufX,fluxBufY,fluxBufZ, lo, &
                                                scaleC, fluxOldX,fluxOldY,fluxOldZ,     &
                                                isFluxDensity)
  use Grid_tile, ONLY : Grid_tile_t
  implicit none

#include "FortranLangFeatures.fh"
  type(Grid_tile_t), intent(in) :: blockDesc
  integer,intent(in) :: lo(3)
  real,intent(in)    :: scaleF,scaleC
  real,intent(in)   ,dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxOldX, fluxOldY, fluxOldZ
  CONTIGUOUS_FSTMT(fluxOldX)
  CONTIGUOUS_FSTMT(fluxOldY)
  CONTIGUOUS_FSTMT(fluxOldZ)
  real,INTENT(INOUT),dimension(size(fluxOldX,1),size(fluxOldX,2),size(fluxOldX,3),size(fluxOldX,4)),TARGET :: fluxBufX
  real,INTENT(INOUT),dimension(size(fluxOldY,1),size(fluxOldY,2),size(fluxOldY,3),size(fluxOldY,4)),TARGET :: fluxBufY
  real,INTENT(INOUT),dimension(size(fluxOldZ,1),size(fluxOldZ,2),size(fluxOldZ,3),size(fluxOldZ,4)),TARGET :: fluxBufZ
  logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate

end subroutine Grid_correctFluxData_xtra
