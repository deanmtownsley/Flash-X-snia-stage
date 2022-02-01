!!****if* source/Grid/GridMain/Grid_getMinCellSize
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
!!  Grid_getMinCellSize
!!
!! SYNOPSIS
!!
!!  call Grid_getMinCellSize(real (OUT)  :: minCellSize)
!!               
!!  
!! DESCRIPTION 
!!
!!  Returns the smallest possible cell size in a simulation in any dimension
!!  that does not represent an angle in curvilinear coordinates.
!!
!!
!! ARGUMENTS
!!
!!  minCellSize - returned value
!!
!!***

subroutine Grid_getMinCellSize(minCellSize)

  use Grid_data, ONLY : gr_minCellSize

  implicit none

  real, intent(OUT) :: minCellSize
  minCellSize = gr_minCellSize

end subroutine Grid_getMinCellSize
