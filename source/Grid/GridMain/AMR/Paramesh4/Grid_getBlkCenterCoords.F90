
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
!!****if* source/Grid/GridMain/paramesh/Grid_getBlkCenterCoords
!!
!! NAME
!!  Grid_getBlkCenterCoords
!!
!! SYNOPSIS
!!
!!  Grid_getBlkCenterCoords(type(block_metadata_t)(IN) :: block,
!!                          real(OUT)   :: blockCenter(MDIM))
!!
!!
!! DESCRIPTION 
!!   Gets the coordinates of the center of the block identified by
!!   blockId.  Returns the coordinates in an array blockCenter
!!
!!
!! ARGUMENTS
!!  block - block metadata
!!  blockCenter - returned array of size MDIM holding the blockCenter coords
!!
!! Example
!!   In 2 dimensions, if physical coordinates are ...
!!    
!!     ________________(0.5 1.0)
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |_______________ |
!!  (-0.5, 0.0)
!!
!!  then the values returned in blockCenter are 
!!  blockCenter(IAXIS) = 0.0
!!  blockCenter(JAXIS) = 0.5
!!  blockCenter(KAXIS) = 0.0 since the dimension is not included  
!!
!!***


subroutine Grid_getBlkCenterCoords(block, blockCenter)

  use tree, ONLY : coord
  use Grid_tile, ONLY : Grid_tile_t
  implicit none

#include "constants.h"

  type(Grid_tile_t),intent(in) :: block
  real,dimension(MDIM),intent(out) :: blockCenter
  blockCenter=coord(:,block%id)

end subroutine Grid_getBlkCenterCoords
