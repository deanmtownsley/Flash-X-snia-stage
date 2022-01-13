!!****if* source/Grid/GridMain/paramesh/Grid_restrictAllLevels
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
!!  Grid_restrictAllLevels
!!
!! SYNOPSIS
!! 
!!  Grid_restrictAllLevels()
!!  
!! DESCRIPTION 
!!  Restricts the grid data to all refinement levels. Normally FLASH
!!  only evolves on the leaf blocks, calling this routine makes all
!!  levels have valid data.  This is mostly for visualization purposes to
!!  be able to look at different levels of resolution
!!  
!!  
!!
!!***


subroutine Grid_restrictAllLevels()

  use Timers_interface, ONLY: Timers_start, Timers_stop

  implicit none
  call Timers_start("restrictAll")
  call gr_restrictTree()
  call Timers_stop("restrictAll")
  
end subroutine Grid_restrictAllLevels
