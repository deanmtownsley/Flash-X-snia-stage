!!****f* source/Grid/Grid_getDeltas
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
!!  Grid_getDeltas
!!
!! SYNOPSIS
!!
!!  call Grid_getDeltas(integer(IN) :: lev,
!!                 real(OUT)   :: del(MDIM))
!!  
!! DESCRIPTION 
!!  
!!  Gets the grid spacing dx/dy/dz for a given level on the Grid.
!!  dx is the size of one cell in the x direction of a block, etc.
!!
!!  
!! ARGUMENTS 
!!
!!  lev - refinement level.
!!        This is 1-based, i.e., the root level is numbered 1.
!!  del - array of size MDIM returned holding the dx, dy, and dz values
!!
!!***

subroutine Grid_getDeltas(lev, del)

  implicit none

#include "constants.h"

  integer, intent(IN)   :: lev
  real, dimension(MDIM), intent(out) :: del
  del=0.0
  return
end subroutine Grid_getDeltas

