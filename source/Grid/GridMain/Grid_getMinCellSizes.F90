!!****if* source/Grid/GridMain/Grid_getMinCellSizes
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!! NAME
!!  Grid_getMinCellSizes
!!
!! SYNOPSIS
!!
!!  Grid_getMinCellSizes(real (OUT)  :: minCellSizes(MDIM))
!!               
!!  
!! DESCRIPTION 
!!
!!  Returns the smallest possible cell size in a simulation for all dimensions.
!!
!!
!! ARGUMENTS
!!
!!  minCellSizes - returned array
!!
!!***

#include "constants.h"

subroutine Grid_getMinCellSizes(minCellSizes)
  use Grid_data, ONLY : gr_minCellSizes
  implicit none
  real, dimension(MDIM), intent(OUT) :: minCellSizes
  minCellSizes(1:MDIM) = gr_minCellSizes(1:MDIM)
end subroutine Grid_getMinCellSizes
