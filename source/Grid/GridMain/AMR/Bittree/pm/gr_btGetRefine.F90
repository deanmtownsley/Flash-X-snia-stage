!!****if* source/Grid/GridMain/paramesh/bittree/source/gr_btGetRefine.F90
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
!!  gr_btGetRefine
!!
!! SYNOPSIS
!!
!!  call gr_btGetRefine(lev,ijk,refine)
!!
!! DESCRIPTION
!!
!!  Checks the refine status of block. Like other amr_bittree routines,
!!  hass normal integer/logical intputs and output. 
!!
!! ARGUMENTS
!!  lev: (in) 1-based level of block
!!  ijk: (in) 0-based block coordinate
!!  derefine: (out) whether block is marked for refine
!!
!!***
subroutine gr_btGetRefine(lev, ijk, refine)
  use bittree, only: gr_btGetBitid, bittree_is_parent, bittree_check_refine_bit
  use Driver_interface, only: Driver_abort
  use iso_c_binding, only: c_int,c_bool

  implicit none
  
  integer, intent(in) :: lev
  integer, intent(in) :: ijk(3)
  logical, intent(out) :: refine
  
  logical(c_bool) :: is_par,marked
  integer :: id, lev1, ijk1(3)

!-Get Bittree ID of block
  lev1 = lev
  ijk1 = ijk
  call gr_btGetBitid(lev1, ijk1, id)

!-Check to make sure block was IDed correctly.
  if ((lev /= lev1).OR.any(ijk/=ijk1))     &
     call Driver_abort("Error identifying block in gr_btGetRefine. &
     &Routine can only be called on existing blocks.")

!-Check to see if block is leaf and marked for nodetype change
  call bittree_is_parent(logical(.FALSE.,c_bool),int(id,c_int),is_par)
  call bittree_check_refine_bit(int(id,c_int),marked)

  refine = (.NOT.is_par).AND.marked

  return
end subroutine
