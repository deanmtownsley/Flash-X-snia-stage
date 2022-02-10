!!****if* source/Grid/GridMain/paramesh/bittree/source/gr_btIsParent.F90
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
!!  gr_btIsParent
!!
!! SYNOPSIS
!!
!!  call gr_btIsParent(lev,ijk,is_par,updated)
!!
!! DESCRIPTION
!!
!!  Checks the parentage of block. Like other amr_bittree routines,
!!  hass normal integer/logical intputs and output. 
!!
!! ARGUMENTS
!!  lev: (in) 1-based level
!!  ijk: (in) 0-based block coordinate
!!  is_par: (out) whether block is parent
!!  updated: (in,optional) true=updated tree, false=original (default=false)
!!
!!***
subroutine gr_btIsParent(lev, ijk, is_par, updated)
  use bittree, only: gr_btGetBitid, bittree_is_parent
  use Driver_interface, only: Driver_abort
  use iso_c_binding, only: c_int,c_bool

  implicit none
  
  integer, intent(in) :: lev
  integer, intent(in) :: ijk(3)
  logical, intent(out) :: is_par
  logical, intent(in), optional :: updated
  
  logical(c_bool) :: check_updated, marked
  integer :: id, lev1, ijk1(3)

  if(present(updated)) then
    check_updated = updated
  else
    check_updated = .FALSE.
  end if

!-Get bittree ID of block
  lev1 = lev
  ijk1 = ijk
  if(present(updated)) then
    call gr_btGetBitid(lev1, ijk1, id, updated)
  else
    call gr_btGetBitid(lev1, ijk1, id )
  end if

!-Check to make sure block was IDed correctly
  if ((lev/=lev1).OR.any(ijk/=ijk1))    then
    call Driver_abort("Error identifying block in gr_btIsParent. &
     &Routine can only be called on existing blocks.")
  end if

!-Check whether parent
  call bittree_is_parent(check_updated,int(id,c_int),marked)
  is_par = marked
  
  return
end subroutine
