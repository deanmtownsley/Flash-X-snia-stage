!!****if* source/Grid/GridMain/paramesh/bittree/source/amr_bittree_get_bitid.F90
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
!!  amr_bittree_get_bitid
!!
!! SYNOPSIS
!!
!!  call amr_bittree_get_bitid(lev,ijk,bitid,updated)
!!
!! DESCRIPTION
!!
!!  Gets the Bittree bit id of a block from its position. This value can
!!  then be used for other functions to check parentage, mark for refinement,
!!  etc. This function has normal integer/logical inputs and outputs.
!!  NOTE: lev and ijk are intent(inout) here, compared to many other functions
!!  where they are only intent(in). This is to allow for this function to be
!!  called to check the existence of a block on a certain level.
!!
!! ARGUMENTS
!!  lev: (inout) 1-based level of block being hunted
!!  ijk: (inout) 0-based block coordinate 
!!  bitid: (out) bit id of block
!!  updated: (in,optional) true=updated tree, false=original (default=false)
!!
!!***
subroutine amr_bittree_get_bitid(lev, ijk, bitid, updated)
  use bittree, only: bittree_identify
  use iso_c_binding, only: c_int,c_bool
  
  implicit none
  
  integer, intent(inout) :: lev
  integer, intent(inout) :: ijk(3)
  integer, intent(out) :: bitid
  logical, intent(in), optional :: updated
  
  integer(c_int) :: mort, levc, ijkc(3), id
  logical(c_bool) :: check_updated

  if(present(updated)) then
    check_updated = updated
  else
    check_updated = .FALSE.
  end if

!-Get 0-based level of block, and c_int version of ijk.
  levc = lev - 1
  ijkc = ijk

!-ID block. The mort output is not used here (use amr_identify_block
!-to return morton number)
  call bittree_identify(check_updated, levc, ijkc, mort, id)

!-Convert back to 1-based level for output.
  lev = levc + 1
  ijk = ijkc
  bitid = id

  return
end subroutine
