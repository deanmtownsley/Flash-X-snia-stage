!!****if* source/Grid/GridMain/AMR/gr_estimateError
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
!!  gr_estimateBlkError
!!
!! SYNOPSIS
!!
!!  call gr_estimateBlkError(real(INOUT) :: error,
!!                   integer(IN) :: iref,
!!                   real(IN)    :: refine_filter)
!!
!!  DESCRIPTION
!!
!!  For one block, estimate the error associated with the given variable to
!!  help determine if the block needs refinement or derefinement.
!!
!!  ARGUMENTS
!!
!!    error - the computed error, a scalar value for the current variable
!!            (given by iref) and current block.
!!
!!    blockDesc - describes the block.
!!
!!    iref - index of the refinement variable in data structure "unk"
!!
!!    refine_filter - makes sure that error calculations to determine refinement
!!                    don't diverge numerically
!!
!!  NOTES
!!
!!    In the case of the PARAMESH Grid implementation, this routine is
!!    called from gr_estimateError.
!!
!!  SEE ALSO
!!
!!    Grid_markRefineDerefine
!!    Grid_markRefineDerefineCallback
!!
!!***

subroutine gr_estimateBlkError(error, tileDesc, iref, refine_filter)
  use Grid_tile, ONLY : Grid_tile_t

  implicit none

  integer, intent(IN) :: iref
  type(Grid_tile_t),intent(IN) :: tileDesc
  real, intent(IN) ::  refine_filter
  real,intent(INOUT) :: error

end subroutine gr_estimateBlkError

