!!****if* source/Grid/GridMain/paramesh/bittree/source/gr_btGetLocalBitids.F90
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
!!  gr_btGetLocalBitids
!!
!! SYNOPSIS
!!
!!  call gr_btGetLocalBitids(mype,bitid_list,updated)
!!
!! DESCRIPTION
!!  Get list of bitids of local blocks.
!!
!! ARGUMENTS
!!  mype: (in) proc of desired processor's blocks
!!  nblks: (in) number of blocks on desired processor
!!  bitid_list: (out) list of blocks' bitids
!!  updated: (in,optional) true=updated tree, false=original tree (treated as false if not present)
!!
!!***
subroutine gr_btGetLocalBitids(mype, nblks, bitid_list, updated)
  use bittree, only: bittree_get_bitid_list, localMortUB, bittree_print_2d, old_localMortUB
  use iso_c_binding, only: c_int,c_bool
  
  implicit none
  
  integer, intent(in) :: mype
  integer, intent(in) :: nblks
  integer, intent(out) :: bitid_list(nblks)
  logical, intent(in), optional :: updated
  
  integer(c_int) :: mmin,mmax
  integer(c_int) :: bitid_listc(nblks)
  logical(c_bool) :: check_updated

  if(present(updated)) then
    check_updated = updated
  else
    check_updated = .FALSE.
  end if


!-Locate block on bittree
  if(mype.eq.0) then
    mmin = 0
  else
    if(check_updated .OR. .NOT.allocated(old_localMortUB)) then
      mmin = localMortUB(mype)
    else
      mmin = old_localMortUB(mype)
    end if
  end if
  if(check_updated .OR. .NOT.allocated(old_localMortUB)) then
    mmax = localMortUB(mype+1)
  else
    mmax = old_localMortUB(mype+1)
  end if
   
  call bittree_get_bitid_list(check_updated, mmin, mmax, bitid_listc)
  bitid_list = int(bitid_listc)

  return
end subroutine
