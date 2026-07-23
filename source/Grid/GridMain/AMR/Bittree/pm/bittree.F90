module bittree
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
  implicit none
  
  public :: localMortUB,old_localMortUB
  public :: gr_btDistributedSort

  logical, save :: gr_btDistributedSort

  integer, allocatable, save :: localMortUB(:)
    !!Highest global morton index for blocks on the specified proc
  integer, allocatable, save :: old_localMortUB(:)
    !!Old version of localMortUB, which disagrees during refinement

  !!!!
  ! General Bittree routines
  !!!!
  
  interface
    subroutine gr_btIdentify(procs, lev, ijk, proc, locblk,updated,bitid)
      integer, intent(in) :: procs
      integer, intent(inout) :: lev
      integer, intent(inout) :: ijk(*)
      integer, intent(out) :: proc
      integer, intent(out) :: locblk
      logical, intent(in), optional :: updated
      integer, intent(out), optional :: bitid
    end subroutine
  end interface
  
  interface
    subroutine gr_btLocate(bitid, lev, ijk, updated)
      integer, intent(in) :: bitid
      integer, intent(out) :: lev
      integer, intent(out) :: ijk(*)
      logical, intent(in), optional :: updated
    end subroutine
  end interface

  interface
    subroutine gr_btGetBitid(lev, ijk, bitid, updated)
      integer, intent(inout) :: lev
      integer, intent(inout) :: ijk(*)
      integer, intent(out) :: bitid
      logical, intent(in), optional :: updated
    end subroutine
  end interface

  interface
    subroutine gr_btIsParent(lev, ijk, is_par, updated)
      integer, intent(in) :: lev
      integer, intent(in) :: ijk(*)
      logical, intent(out) :: is_par
      logical, intent(in), optional :: updated
    end subroutine
  end interface

  interface
    subroutine gr_btRefineMark(lev, ijk, val)
      integer, intent(in) :: lev
      integer, intent(in) :: ijk(*)
      logical, intent(in), optional :: val
    end subroutine
  end interface
  
  interface
    subroutine gr_btDerefineMark(lev, ijk, val)
      integer, intent(in) :: lev
      integer, intent(in) :: ijk(*)
      logical, intent(in), optional :: val
    end subroutine
  end interface

  interface
    subroutine gr_btGetRefine(lev, ijk, refine)
      integer, intent(in) :: lev
      integer, intent(in) :: ijk(*)
      logical, intent(out) :: refine
    end subroutine
  end interface
  
  interface
    subroutine gr_btGetDerefine(lev, ijk, derefine)
      integer, intent(in) :: lev
      integer, intent(in) :: ijk(*)
      logical, intent(out) :: derefine
    end subroutine
  end interface

  interface
    subroutine gr_btGetLocalBitids(mype, nblks, bitid_list,updated)
      integer, intent(in) :: mype
      integer, intent(in) :: nblks
      integer, intent(out) :: bitid_list(*)
      logical, intent(in), optional :: updated
    end subroutine
  end interface


  !!!!!!!
  ! These need to be implemented by the Amr package, currently only done for PM
  !!!!!!!

  interface
    subroutine gr_btSortMortonBittree(nprocs,mype,sort_by_work)
      integer, intent(in) :: nprocs
      integer, intent(in) :: mype
      logical, intent(in),optional :: sort_by_work
    end subroutine
  end interface

  interface
    subroutine gr_getIntCoords(lblock,lcoord)
      integer, intent(in) :: lblock
      integer, intent(out), dimension(3) :: lcoord
    end subroutine
  end interface

  interface
    subroutine gr_getNeighIntCoords(lblock,gCell,neighCoord,force)
      integer, intent(in) :: lblock
      integer, intent(in), dimension(3) :: gCell
      integer, intent(out), dimension(3) :: neighCoord
      logical, intent(in), optional :: force
    end subroutine
  end interface
  
  interface
    subroutine amr_build_bittree
    end subroutine amr_build_bittree
  end interface

  interface
    subroutine amr_verify_bittree
    end subroutine
  end interface

  interface
    subroutine amr_calculate_tree_data(nprocs,mype,lnblocks_old)
     integer, intent(in) :: nprocs
     integer, intent(in) :: mype
     integer, intent(in) :: lnblocks_old
    end subroutine
  end interface

  interface
    subroutine amr_morton_order_bittree(nprocs,mype,ref_count)
      integer, intent(in) :: nprocs
      integer, intent(in) :: mype
      integer, intent(in) :: ref_count
    end subroutine
  end interface

  interface
    subroutine amr_exchange_work_bflags(nprocs,mype,lnblocks_old,new_loc,old_loc,new_child)
      integer, intent(in) :: nprocs,mype
      integer, intent(in) :: lnblocks_old
      integer, intent(in)  :: new_loc(*)
      integer, intent(in)  :: old_loc(*)
      integer, intent(in)  :: new_child(*)
    end subroutine
  end interface

  !!!!!
  ! Fortran C Interfaces
  !!!!!
 
  interface
    function bittree_initialized() result(yep) bind(c,name='bittree_initialized')
      use iso_c_binding, only: c_bool
      logical(c_bool) :: yep
    end function
  end interface
  
  interface
    subroutine bittree_init(ndim, top, topmask) bind(c,name='bittree_init')
      use iso_c_binding, only: c_bool, c_int
      integer(c_int), intent(in) :: ndim
      integer(c_int), intent(in) :: top(*)
      logical(c_bool), intent(in) :: topmask(*)
    end subroutine
  end interface
  
  interface
    subroutine bittree_block_count(updated, bcount) bind(c,name='bittree_block_count')
      use iso_c_binding, only: c_bool, c_int
      logical(c_bool),intent(in)  :: updated
      integer(c_int), intent(out) :: bcount
    end subroutine
  end interface
  
  interface
    subroutine bittree_leaf_count(updated, lcount) bind(c,name='bittree_leaf_count')
      use iso_c_binding, only: c_bool, c_int
      logical(c_bool),intent(in)  :: updated
      integer(c_int), intent(out) :: lcount
    end subroutine
  end interface
  
  interface
    subroutine bittree_delta_count(dcount) bind(c,name='bittree_delta_count')
      use iso_c_binding, only: c_bool, c_int
      integer(c_int), intent(out) :: dcount
    end subroutine
  end interface

  interface
    subroutine bittree_check_refine_bit(bitid, bitval) bind(c,name='bittree_check_refine_bit')
      use iso_c_binding, only: c_bool,c_int
      integer(c_int), intent(in) :: bitid
      logical(c_bool),intent(out) :: bitval
    end subroutine
  end interface
 
  interface
    subroutine bittree_is_parent(updated, bitid, parent) bind(c,name='bittree_is_parent')
      use iso_c_binding, only: c_int,c_bool
      logical(c_bool), intent(in) :: updated
      integer(c_int), intent(in) :: bitid
      logical(c_bool), intent(out) :: parent
    end subroutine
  end interface
  
  interface
    subroutine bittree_identify(updated, lev, ijk, mort, bitid) bind(c,name='bittree_identify')
      use iso_c_binding, only: c_int,c_bool
      logical(c_bool), intent(in)   :: updated
      integer(c_int), intent(inout) :: lev
      integer(c_int), intent(inout) :: ijk(*)
      integer(c_int), intent(out) :: mort
      integer(c_int), intent(out) :: bitid
    end subroutine
  end interface
  
  interface
    subroutine bittree_locate(updated, bitid, lev, ijk, mort) bind(c,name='bittree_locate')
      use iso_c_binding, only: c_int,c_bool
      logical(c_bool), intent(in)   :: updated
      integer(c_int), intent(in) :: bitid
      integer(c_int), intent(out) :: lev
      integer(c_int), intent(out) :: ijk(*)
      integer(c_int), intent(out) :: mort
    end subroutine
  end interface

  interface
    subroutine bittree_get_bitid_list(updated, mort_min, mort_max, idout) bind(c,name='bittree_get_bitid_list')
      use iso_c_binding, only: c_int,c_bool
      logical(c_bool), intent(in)   :: updated
      integer(c_int), intent(in) :: mort_min
      integer(c_int), intent(in) :: mort_max
      integer(c_int), intent(out) :: idout(*)
    end subroutine
  end interface

  interface
    subroutine bittree_refine_init() bind(c,name='bittree_refine_init')
    end subroutine
  end interface

  interface
    subroutine bittree_refine_mark(bitid, val) bind(c,name='bittree_refine_mark')
      use iso_c_binding, only: c_int,c_bool
      integer(c_int), intent(in) :: bitid
      logical(c_bool), intent(in) :: val
    end subroutine
  end interface

  interface
    subroutine bittree_refine_reduce(comm) bind(c,name='bittree_refine_reduce')
      use iso_c_binding, only: c_int
      integer(c_int), intent(in) :: comm
    end subroutine
  end interface
  
  interface
    subroutine bittree_refine_reduce_and(comm) bind(c,name='bittree_refine_reduce_and')
      use iso_c_binding, only: c_int
      integer(c_int), intent(in) :: comm
    end subroutine
  end interface
  
  interface
    subroutine bittree_refine_update() bind(c,name='bittree_refine_update')
    end subroutine
  end interface
  
  interface
    subroutine bittree_refine_apply() bind(c,name='bittree_refine_apply')
    end subroutine
  end interface
  
  interface
    subroutine bittree_print_2d(datatype) bind(c,name='bittree_print_2d')
      use iso_c_binding, only: c_int
      integer(c_int), intent(in) :: datatype
    end subroutine
  end interface
end module
