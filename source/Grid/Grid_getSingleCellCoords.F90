!!****f* source/Grid/Grid_getSingleCellCoords
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
!!  Grid_getSingleCellCoords
!!
!! SYNOPSIS
!!
!!  call Grid_getSingleCellCoords(integer(in)       ::  ind(MDIM), 
!!                           integer(in)       ::  level,
!!                           integer(in)       ::  edge,
!!                           real(OUT)         ::  coords(MDIM))
!!  
!! DESCRIPTION 
!!
!!  Returns the coordinates of a single cell of a given block.
!!
!!
!! ARGUMENTS
!!
!!  ind - array holding the indices of the cell whose coordinates to return.
!!        In this implementation, global indices (for the given level) are assumed.
!! 
!!  level - refinement level (1 based)
!!
!!  edge - indicates if user wants the left, center or right edge of the block
!!         options are LEFT_EDGE, RIGHT_EDGE, or CENTER.  These constants are
!!         defined in constants.h
!!
!!    
!!  coords : returned coordinates of the specificed cell
!!
!!***

#include "constants.h"

subroutine Grid_getSingleCellCoords(ind, level, edge, coords)
  use Driver_interface, ONLY : Driver_abort

  implicit none

  integer, intent(in)  :: ind(1:MDIM)
  integer, intent(in)  :: level
  integer, intent(in)  :: edge
  real,    intent(out) :: coords(1:MDIM)

  coords(:) = 0.0

  call Driver_abort("[Grid_getSingleCellCoords] DEPRECATED")
end subroutine Grid_getSingleCellCoords
