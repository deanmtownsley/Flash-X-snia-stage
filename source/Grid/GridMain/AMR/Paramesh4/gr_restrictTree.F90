!!****if* source/Grid/GridMain/paramesh/gr_restrictTree
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
!!  gr_restrictTree
!!
!! SYNOPSIS
!!
!!  gr_restrictTree()
!!
!!  
!! DESCRIPTION 
!!
!!  This subroutine restrict the tree all the way to the coarsest
!!  refinement level. It is usually invoked before writing data to
!!  and output file. This is mostly for visualization purposes to
!!  be able to look at different levels of resolution
!!
!!
!! ARGUMENTS   
!!
!!***

subroutine gr_restrictTree()

  use physicaldata, ONLY : no_permanent_guardcells
  use tree, ONLY : nodetype,newchild,lnblocks,lrefine,maxblocks_tr

  use paramesh_interfaces, only : amr_restrict
  
  use Grid_data, ONLY:  gr_meshComm, gr_meshMe
  
#include "Flashx_mpi_implicitNone.fh"
#include "constants.h"
  integer :: gridDataStruct
  
  integer nodetype_old(maxblocks_tr),level_at,lmax,lmax2
  integer i,ierr
  logical newchild_old(maxblocks_tr)

  lmax = -1
  do i = 1,lnblocks
     if (nodetype(i).eq.1) then
        lmax = max(lrefine(i),lmax)
     end if
  end do
  
  call MPI_ALLREDUCE (lmax,lmax2,1,MPI_INTEGER, &
       &     MPI_MAX,gr_meshComm,ierr)
  lmax = lmax2
  
  nodetype_old(:) = nodetype(:)
  newchild_old(:) = newchild(:)
  
  do level_at = lmax,2,-1
     
     call amr_restrict(gr_meshMe,1,0,.false.)
     
     do i = 1, lnblocks
        if (nodetype(i).eq.1.and.lrefine(i).eq.level_at) &
             &           nodetype(i) = -1
        if (nodetype(i).eq.2.and.lrefine(i).eq.level_at-1) &
             &           nodetype(i) = 1
        if (nodetype(i).eq.3.and.lrefine(i).eq.level_at-2) &
             &           nodetype(i) = 2
     end do
     
     do i = 1,lnblocks
        if (nodetype(i).eq.100) nodetype(i) = 1
     end do
     
     ! reset neighbor nodetypes and child nodetypes
     call amr_morton_process()
     
  end do
  
  nodetype(:) = nodetype_old(:)
  newchild(:) = newchild_old(:)
  ! reset neighbor nodetypes and child nodetypes
  
  call amr_morton_process()

  if (no_permanent_guardcells) then
     gridDataStruct=CENTER
#if NFACE_VARS > 0
     gridDataStruct=CENTER_FACES
#endif
     call gr_commSetup(gridDataStruct)
  else
     call gr_freeCommRecvBuffer
  end if

  return
end subroutine gr_restrictTree
