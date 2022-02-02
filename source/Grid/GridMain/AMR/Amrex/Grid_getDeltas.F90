!!****if* source/Grid/GridMain/AMR/Amrex/Grid_getDeltas
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
!!  Grid_getDeltas(integer(IN) :: level,
!!                 real(OUT)   :: del(MDIM))
!!  
!! DESCRIPTION 
!!  Gets the size (dX, dY, dZ) of cells in the given refinement level where dX
!!  is the width of the cell along the first dimension (IAXIS).
!! 
!!  If a dimension is Cartesian, then the returned size is the length of any
!!  cell along that direction.  If a dimension is angular, then the size is the
!!  angle expressed in radians (as opposed to the arclength).
!!
!! ARGUMENTS 
!!  level - local block number
!!  del - array of size MDIM returned holding the dX, dY, and dZ values
!!
!!***

subroutine Grid_getDeltas(lev, del)
  use amrex_amrcore_module, ONLY : amrex_geom
  
  implicit none

#include "constants.h"
  
  integer, intent(IN)  :: lev
  real,    intent(OUT) :: del(MDIM)

  ! AMReX uses zero-based level indexing, but FLASH assumes one-based
  del = 0.0
  del(1:MDIM) = amrex_geom(lev-1)%dx(1:MDIM)
end subroutine Grid_getDeltas

