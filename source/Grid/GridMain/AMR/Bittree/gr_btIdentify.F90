!!****if* source/Grid/GridMain/paramesh/bittree/source/gr_btIdentify
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
!!  gr_btIdentify
!!
!! SYNOPSIS
!!
!!  call gr_btIdentify(procs,lev,ijk,proc,locblk,updated,bitid)
!!
!! DESCRIPTION
!!
!!  Determines the processor a block resides on from its spatial index.
!!  If a block matching the lev,ijk coordinates exists then its proc,locblk are
!!  filled in accordingly and lev,ijk are returned untouched.  If that block
!!  doesnt exist, then the finest parent containing that space will be returned
!!  in lev,ijk,proc,locblk. Hence block existence can be determined from change
!!  in the lev argument. If that space lies outside the domain then lev=0
!!  will be returned, boundary periodicity is not considered here. To consider 
!!  boundary periodicity, use gr_getNeighIntCoords to get neighbor coordinates.
!!
!! ARGUMENTS
!!  procs: (in) number of processors on mesh (result of MPI_COMM_SIZE)
!!  lev: (inout) 1-based level of block being hunted
!!  ijk: (inout) 0-based block coordinate (see gr_xyzToBlockLevel)
!!  proc: (out) 0-based processor rank block resides on
!!  locblk: (out) 1-based local block index for that processor
!!  updated: (in,optional) true=updated tree, false=original tree (treated as false if not present)
!!  bitid: (out,optional) output of the bitid of desired bloc 
!!
!!***
subroutine gr_btIdentify(procs, lev, ijk, proc, locblk, updated, bitid)
  use gr_sortByWorkTools, only: gr_calcProcLocblk
  use bittree, only: bittree_block_count, bittree_identify, localMortUB,old_localMortUB
  use iso_c_binding, only: c_int,c_bool
  
  implicit none
  
  integer, intent(in) :: procs
  integer, intent(inout) :: lev
  integer, intent(inout) :: ijk(3)
  integer, intent(out) :: proc
  integer, intent(out) :: locblk
  logical, intent(in), optional :: updated
  integer, intent(out), optional :: bitid
  
  integer(c_int) :: levc, ijkc(3), id, blks, mort
  logical(c_bool) :: check_updated
  integer :: j, mort1

  if(present(updated)) then
    check_updated = updated
  else
    check_updated = .FALSE.
  end if


!-Identify block on bittree
  levc = lev - 1
  ijkc = ijk
  call bittree_identify(check_updated, levc, ijkc, mort, id)
  lev = levc + 1
  ijk = ijkc
  if(present(bitid)) bitid = int(id)

  mort1 = int(mort) + 1

!-Calculate proc and locblk.
  if (check_updated .OR. .NOT.allocated(old_localMortUB)) then
    call gr_calcProcLocblk(procs,localMortUB,mort1,proc,locblk)
  else
    call gr_calcProcLocblk(procs,old_localMortUB,mort1,proc,locblk)
  end if
  
  return
end subroutine
