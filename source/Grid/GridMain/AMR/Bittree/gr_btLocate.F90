!!****if* source/Grid/GridMain/paramesh/bittree/source/gr_btLocate.F90
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
!!  gr_btLocate
!!
!! SYNOPSIS
!!
!!  call gr_btLocate(bitid,lev,ijk,updated)
!!
!! DESCRIPTION
!!  Locate a block (get its lrefine and integer coordinates) from its bittree id.
!!
!! ARGUMENTS
!!  bitid: (in) the bitid of desired block 
!!  lev: (out) 1-based level of block being hunted
!!  ijk: (out) 0-based block coordinate 
!!  updated: (in,optional) true=updated tree, false=original tree (treated as false if not present)
!!
!!***
#include "Simulation.h"
subroutine gr_btLocate(bitid, lev, ijk, updated)
  use bittree, only: bittree_block_count, bittree_locate
  use iso_c_binding, only: c_int,c_bool
  
  implicit none
  
  integer, intent(in) :: bitid
  integer, intent(out) :: lev
  integer, intent(out) :: ijk(3)
  logical, intent(in), optional :: updated
  
  integer(c_int) :: levc, ijkc(3), mort
  logical(c_bool) :: check_updated

  if(present(updated)) then
    check_updated = updated
  else
    check_updated = .FALSE.
  end if


!-Locate block on bittree 
  ijk(:) = 0
  call bittree_locate(check_updated, int(bitid,c_int), levc, ijkc, mort)
  lev = levc + 1
  ijk(1:NDIM) = ijkc(1:NDIM)

  return
end subroutine
