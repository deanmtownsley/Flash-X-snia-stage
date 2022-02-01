!!****f* source/Grid/Grid_restrictAllLevels
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
!!  Grid_restrictAllLevels
!!
!! SYNOPSIS
!! 
!!  call Grid_restrictAllLevels()
!!  
!! DESCRIPTION 
!!  Restricts the grid data to all refinement levels. Normally FLASH
!!  only evolves on the leaf blocks, calling this routine makes all
!!  levels have valid data.  This is mostly for visualization purposes to
!!  be able to look at different levels of resolution
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  For the Uniform Grid or any other mesh package with only a single 
!!  level, no source file at the implementation level of the Grid unit
!!  structure is provided, so that the stub implementation will be used.
!!
!!***


subroutine Grid_restrictAllLevels()

implicit none
end subroutine Grid_restrictAllLevels
