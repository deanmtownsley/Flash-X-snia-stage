!!****if* source/Grid/GridMain/AMR/Amrex/gr_estimateError
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
!!  gr_estimateError
!!
!! SYNOPSIS
!!
!!  call gr_estimateError(real(INOUT) :: error(MAXBLOCKS),
!!                   integer(IN) :: iref,
!!                   real(IN)    :: refine_filter)
!!  
!!  DESCRIPTION
!!  
!!  For each block, estimate the error associated with the given variable to
!!  help determine if the block needs refinement or derefinement.  Update the
!!  corresponding value in the error array to be the maximum of the incoming
!!  value and the value calculated here.
!!
!!  ARGUMENTS 
!!
!!    error - indexed by block IDs.
!!
!!    iref - index of the refinement variable in data structure "unk"
!!
!!    refine_filter - makes sure that error calculations to determine refinement
!!                    don't diverge numerically 
!! 
!!  NOTES
!!  
!!    See Grid_markRefineDerefine
!!
!!    DEV: Not usable with the Amrex Grid implementation, will abort;
!!    however, gr_estimateBlkError is available with the Amrex Grid.
!!
!!  SEE ALSO
!!  
!!    Grid_markRefineDerefine
!!    gr_estimateBlkError
!!
!!***

subroutine gr_estimateError(error, iref, refine_filter)
  use Driver_interface, ONLY : Driver_abort

  implicit none

  integer, intent(IN)    :: iref
  real,    intent(IN)    :: refine_filter
  real,    intent(INOUT) :: error(MAXBLOCKS)
 
  call Driver_abort("[gr_estimateError] Not implemented with AMReX")
end subroutine gr_estimateError

