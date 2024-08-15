!!****if* source/Grid/GridMain/paramesh/bittree/source/gr_btGetDerefine.F90
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
!!  gr_btGetDerefine
!!
!! SYNOPSIS
!!
!!  call gr_btGetDerefine(lev,ijk,derefine)
!!
!! DESCRIPTION
!!
!!  Checks the derefine status of block. Like other amr_bittree routines,
!!  hass normal integer/logical intputs and output. 
!!
!! ARGUMENTS
!!  lev: (in) 1-based level of block
!!  ijk: (in) 0-based block coordinate
!!  derefine: (out) whether block is marked for derefine
!!
!!***
subroutine gr_btGetDerefine(lev, ijk, derefine)
  use bittree, only: gr_btGetBitid, bittree_is_parent, bittree_check_refine_bit
  use Driver_interface, only: Driver_abort
  use iso_c_binding, only: c_int,c_bool
  implicit none
  
  integer, intent(in) :: lev
  integer, intent(in) :: ijk(3)
  logical, intent(out) :: derefine
  
  logical(c_bool) :: is_par,marked
  integer :: id, lev_par, ijk_par(3)

  if(lev.eq.1) then
    derefine=.FALSE.
    return
  endif

!-Get lev and ijk of parent
  lev_par = lev - 1
  ijk_par = ijk/2
 
!-Get bittree ID of parent
  call gr_btGetBitid(lev_par, ijk_par, id)

!-Check to make sure parent was IDed correctly
  if ((lev /= lev_par+1).OR.any(ijk/=(ijk_par*2+mod(ijk,2))))  &
     call Driver_abort("Error identifying block in gr_btGetDerefine. &
     &Routine can only be called on existing blocks.")

!-Make sure parent is both marked as a parent and for nodetype change
  call bittree_is_parent(logical(.FALSE.,c_bool), int(id,c_int), is_par)
  call bittree_check_refine_bit(int(id,c_int), marked)

  derefine = (is_par).AND.marked
  return
end subroutine
