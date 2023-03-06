!!****if* source/Grid/localAPI/gr_markRefineDerefine
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  gr_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  call gr_markRefineDerefine(real(IN) :: error(MAXBLOCKS),
!!                             real(IN) :: refine_cutoff,
!!                             real(IN) :: derefine_cutoff)
!!  
!!  DESCRIPTION
!!  
!!    Blocks are marked for refining or derefining based on error estimates.
!!
!!    The error estimates passed in as the first argument are typically obtained by
!!    second derivative calculations on one of a set of specified variables to
!!    determine if the block needs more resolution (refine) or less resolution (derefine).
!!    The arguments de/refine_cutoff are the thresholds for triggering the corresponding action.
!!
!!    After blocks have been marked, control is meant to be passed to Paramesh for actually
!!    updating the refinement of the grid.
!!
!!  ARGUMENTS
!!
!!    error - an array containing estimates of the largest errors of each block.
!!            These numbers are relative numbers typically in the range 0.0 .. 1.0
!!            that are directly comparable to refine_cutoff and derefine_cutoff values.
!!
!!    refine_cutoff - the threshold value for triggering refinement 
!!
!!    derefine_cutoff - the threshold for triggering derefinement
!!
!!  NOTES
!!  
!!    See Grid_markRefineDerefine
!!
!!  SEE ALSO
!!
!!    gr_estimateError
!!    Grid_markRefineDerefine
!!
!!***

subroutine gr_markRefineDerefine(error, refine_cutoff, derefine_cutoff)

  implicit none

#include "Simulation.h"

  real, intent(IN) :: error(MAXBLOCKS)
  real, intent(IN) :: refine_cutoff, derefine_cutoff

end subroutine gr_markRefineDerefine

