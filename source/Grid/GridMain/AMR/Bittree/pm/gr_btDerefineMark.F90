!!****if* source/Grid/GridMain/paramesh/bittree/source/gr_btDerefineMark.F90
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
!!
!!  gr_btDerefineMark
!!
!! SYNOPSIS
!!
!!  call gr_btDerefineMark(lev,ijk,val)
!!
!! DESCRIPTION
!!
!!  Mark or unmarks parent for nodetype change in Bittree. In effect, this
!!  marks block at lev, ijk for derefinement if val=True. Like other
!!  amr_bittree routines, this function take regular integer/logical inputs,
!!  and converts them to c_int/c_bool for passing to C interface functions.
!!
!! ARGUMENTS
!!  lev: (in) 1-based level of block 
!!  ijk: (in) 0-based block coordinate 
!!  val: (in,optional) value to set block in delta tree, default = True
!!
!!***
subroutine gr_btDerefineMark(lev, ijk, val)
  use bittree, only: gr_btGetBitid, bittree_refine_mark,bittree_is_parent
  use iso_c_binding, only: c_int,c_bool
  use Driver_interface, ONLY : Driver_abort

  implicit none
  
  integer, intent(in) :: lev
  integer, intent(in) :: ijk(3)
  logical, intent(in),optional :: val
  
  integer :: lev_par, ijk_par(3), id
  logical(c_bool) :: mark, is_par

  if(present(val)) then
    mark = val
  else
    mark = .TRUE.
  end if

!-Compute lrefine and coordinates of parent
  lev_par = lev - 1
  ijk_par = ijk/2

!-Get bittree ID of parent
  call gr_btGetBitid(lev_par, ijk_par, id)
  
!-Make sure bittree identified correct block
  if ((lev /= lev_par + 1).OR.any(ijk /= ijk_par*2 + mod(ijk,2))) &
   call Driver_abort("Error identifying block in gr_btDerefineMark. &
     &Routine can only be called on existing blocks.")

!-Abort if trying to mark leaf
  call bittree_is_parent(logical(.FALSE.,c_bool),int(id,c_int),is_par)
  if (.NOT.is_par) &
    call Driver_abort("Error in gr_btDerefineMark. &
    &Trying to mark leaf block for derefinement (bittree marks parents for derefinement).")
 
!-Mark parent on Bittree's refine delta
  call bittree_refine_mark(int(id,c_int), mark)

  return
end subroutine
