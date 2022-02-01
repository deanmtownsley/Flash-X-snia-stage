!!****if* source/Grid/GridMain/AMR/Grid_markRefineDerefine
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
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  A stub implementation as this public Grid routine is not needed for
!!  AMReX-based FLASH simulations.
!!
!!  Refinement of blocks is managed by AMReX through through the
!!  gr_markRefineDerefineCallback routine.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine Grid_markRefineDerefine()
  use Driver_interface, ONLY : Driver_abort
  
  implicit none
  
  call Driver_abort("[Grid_markRefineDerefine] Not implemented in AMReX")
end subroutine Grid_markRefineDerefine

