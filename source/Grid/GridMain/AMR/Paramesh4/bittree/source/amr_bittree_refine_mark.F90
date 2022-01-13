!!****if* source/Grid/GridMain/paramesh/bittree/source/amr_bittree_refine_mark.F90
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
!!
!!  amr_bittree_refine_mark

!!
!! SYNOPSIS
!!
!!  call amr_bittree_refine_mark(lev,ijk,val)
!!
!! DESCRIPTION
!!
!!  Mark or unmarks leafe block for nodetype change in Bittree. In effect, this
!!  marks block at lev, ijk for refinement if val=True. Like other
!!  amr_bittree routines, this function take regular integer/logical inputs,
!!  and converts them to c_int/c_bool for passing to C interface functions.
!!
!! ARGUMENTS
!!  lev: (in) 1-based level of block
!!  ijk: (in) 0-based block coordinate
!!  val: (in,optional) value to set block in delta tree, default = True
!!
!!***
subroutine amr_bittree_refine_mark(lev, ijk, val)
  use bittree, only: amr_bittree_get_bitid, bittree_refine_mark, bittree_is_parent
  use Driver_interface, only: Driver_abortFlash
  use iso_c_binding, only: c_int,c_bool

  implicit none
  
  integer, intent(in) :: lev
  integer, intent(in) :: ijk(3)
  logical, intent(in), optional :: val
  
  logical(c_bool) :: mark, is_par
  integer :: id, lev1, ijk1(3)

  if (present(val)) then
    mark = val
  else
    mark = .TRUE.
  end if

!-Get bittree ID of block.
  lev1 = lev
  ijk1 = ijk
  call amr_bittree_get_bitid(lev1, ijk1, id)

!-Check to make sure block was IDed correctly.
  if ((lev /= lev1).OR.any(ijk /= ijk1)) &
   call Driver_abortFlash("Error identifying block in amr_bittree_refine_mark. &
     &Routine can only be called on existing blocks.")

!-Abort if trying to mark parent.
  call bittree_is_parent(logical(.FALSE.,c_bool),int(id,c_int),is_par)
  if (is_par) call Driver_abortFlash("Error in amr_bittree_refine_mark. &
    &Trying to mark parent block for refinement")

!-Updated block on bittree's refine delta.
  call bittree_refine_mark(int(id,c_int), mark)

end subroutine
