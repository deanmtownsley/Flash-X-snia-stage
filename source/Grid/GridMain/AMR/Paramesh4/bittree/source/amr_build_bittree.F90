!> \file
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @brief call amr_build_bittree()
!! NAME
!!
!!  amr_build_bittree
!!
!! SYNOPSIS
!!
!!  call amr_build_bittree()
!!
!! DESCRIPTION
!!
!!  Constructs the bittree from scratch by scanning all blocks in the mesh
!!  communicator.
!!
#include "constants.h"

subroutine amr_build_bittree()
  use bittree, only: bittree_init,amr_bittree_refine_mark, &
                  bittree_block_count, bittree_refine_init,&
                  bittree_refine_update,bittree_refine_reduce,&
                  bittree_refine_apply,amr_verify_bittree,&
                  amr_sort_morton_bittree, gr_getIntCoords, &
                  gr_btDistributedSort, gr_btExchangeBflags
  use paramesh_dimensions, only: ndim
  use Paramesh_comm_data, only: amr_mpi_meshComm
  use tree, only: lnblocks, coord, bsize, lrefine, nodetype, &
                  grid_xmin, grid_ymin, grid_zmin, &
                  grid_xmax, grid_ymax, grid_zmax, &
                  gr_btSortByWork,gr_btCustomWork,&
                  gr_btWorkDefaultLeaf,gr_btWorkDefaultPar, &
                  gr_btWorkBoundsPar, gr_btWorkBoundsLeaf

  use Driver_interface, only: Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_interface, only : Grid_setWorkBounds
  use gr_specificData, ONLY : gr_gidIsValid

  use iso_c_binding, only: c_bool, c_int
  
  implicit none
#include "mpif.h"
  
  real :: gmin(ndim), gmax(ndim)
  integer :: top(3)
  integer :: ijk(3),lev
  logical, allocatable :: topmask(:,:,:)
  integer(c_int) :: b, pop0, pop1
  integer :: ierr, nprocs, mype
  logical (c_bool) :: cfalse = .FALSE.

  
  Call MPI_COMM_SIZE (amr_mpi_meshComm,nprocs,ierr) 
  Call MPI_COMM_RANK (amr_mpi_meshComm,mype,ierr) 

!-----------------------------------------------------------------------------------------------
!-Create top(3), where top(1)=Nblockx, etc.
  gmin = reshape((/grid_xmin, grid_ymin, grid_zmin/), (/ndim/))
  gmax = reshape((/grid_xmax, grid_ymax, grid_zmax/), (/ndim/))
  !all procs with blocks *should* agree on the nblockx/y/z values
  if(lnblocks > 0) then
    top(:) = 1
    top(1:ndim) = ishft(int(floor(0.5 + (gmax-gmin)/bsize(1:ndim,1))), 1-lrefine(1))
  else
    top(:) = 0
  end if
  call MPI_ALLREDUCE(MPI_IN_PLACE, top, 3, MPI_INTEGER, MPI_BOR, amr_mpi_meshComm, ierr)

!-----------------------------------------------------------------------------------------------
!-Create topmask, which marks which root blocks exist. For now, if any do not exist, abort.
  allocate(topmask(top(1),top(2),top(3)))
  topmask(:,:,:) = .false.
  do b=1, lnblocks
    if(lrefine(b) == 1) then
      call gr_getIntCoords(b,ijk)
      ijk = 1+ijk
      !ijk(1:ndim) = 1 + int((coord(1:ndim,b) - gmin)/bsize(1:ndim,b))
      topmask(ijk(1),ijk(2),ijk(3)) = .true.
    end if
  end do
  call MPI_ALLREDUCE(MPI_IN_PLACE, topmask, product(top), MPI_LOGICAL, MPI_LOR, amr_mpi_meshComm, ierr)
  
  if(.NOT.all(topmask))     &
    call Driver_abortFlash("Error in amr_build_bittree. All possible root blocks must exist.")

!-Set some runtime parameters needed for amr_sort_morton_bittree
  call RuntimeParameters_get("gr_btDistributedSort",gr_btDistributedSort)
  call RuntimeParameters_get("gr_btExchangeBflags",gr_btExchangeBflags)
  if(gr_btDistributedSort.AND..NOT.(gr_btSortByWork.AND.gr_btCustomWork)) &
    print *,"Warning: For bittree runtime parameters. gr_btDistributedSort &
             &only has effect if both gr_btSortByWork and gr_btCustomWork are True"

!-----------------------------------------------------------------------------------------------
!-create bittree
  call bittree_init(ndim, int(top,c_int), logical(topmask,c_bool))
  
!-----------------------------------------------------------------------------------------------
!-Create finer levels than top on bittree.
  lev = 1
  pop0 = 0
  call bittree_block_count(cfalse, pop1)
!-Loop as long as the previous loop created new blocks on bittree
  do while(pop0 < pop1)
    call bittree_refine_init()
!---If blocks have nodetype>1, mark for refinement
    do b=1, lnblocks
      if(lrefine(b) == lev .and. nodetype(b) > 1) then
        call gr_getIntCoords(b,ijk)
        call amr_bittree_refine_mark(lev, ijk)
      end if
    end do
!---Apply refinement to bittree
    call bittree_refine_reduce(amr_mpi_meshComm)
    call bittree_refine_update()
    call bittree_refine_apply()
    pop0 = pop1
    call bittree_block_count(cfalse,pop1)
    lev = lev + 1
  end do

  call amr_sort_morton_bittree(nprocs,mype,sort_by_work=.FALSE.)
  if (.NOT. gr_gidIsValid) then
    ! bittree check fails as restart file is a AMReX file
  else
    call amr_verify_bittree()
  endif
  ! 
end subroutine
