!!****if* source/Grid/GridMain/paramesh/Grid_receiveInputData
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
!!  Grid_receiveInputData
!!
!! SYNOPSIS
!!
!!  call Grid_receiveInputData(integer(IN) :: localNumBlocks,
!!                             integer(IN) :: alnblocks,
!!                             integer(IN) :: xx)
!!
!! DESCRIPTION 
!!
!!  Initializes grid arrays from arrays read by the I/O unit.
!!
!! ARGUMENTS  
!!
!!  localNumBlocks : the number of blocks on my processor.
!!
!!  alnblocks : the approximate number of local blocks on each
!!              processor if we give each processor an equal
!!              number of blocks.  Calculated from
!!              int(globalNumBlocks/meshNumProcs) + 1.
!!
!!  xx : an integer representing a cutoff point.  Processors
!!       less than this value are assigned alnblocks blocks and
!!       processors greater than or equal to this value are
!!       assigned lnblocks-1 blocks.
!!
!! NOTES
!!
!!  The Paramesh4 implementation initializes the following arrays
!!  owned by PARAMESH:
!!     neigh
!!     parent
!!     child
!!     surr_blks     (conditionally)
!!  using data from the following arrays (presumably filled by io_readData):
!!     gr_gid
!!     gr_gsurr_blks (conditionally)
!!
!!  Other PARAMESH-owned data will have to be initialized elsewhere.
!!
!!  Grid implementations other than Paramesh4 will have to perform
!!  initialization of their internal data structures elswhere.
!!***

subroutine Grid_receiveInputData(localNumBlocks, alnblocks, xx)

#include "constants.h"
#include "Simulation.h"

  use tree, ONLY : nfaces, nchild, neigh, parent, child
#ifdef FLASH_GRID_PARAMESH3OR4
  use gr_specificData, ONLY : gr_gsurr_blks, gr_is_gsurr_blks_initialized
  use physicaldata, ONLY : gsurrblks_set
#endif
#ifdef FLASH_GRID_PARAMESH4DEV_SURR_BLKS_OPTION
  use tree, ONLY : surr_blks
  use physicaldata, ONLY : use_flash_surr_blks_fill, surr_blks_valid
#endif
  use gr_specificData, ONLY : gr_gid, gr_gidIsValid
  use Grid_Data, ONLY:  gr_globalMe, gr_meshNumProcs

  implicit none
  integer, intent(IN) :: localNumBlocks, alnblocks, xx
  integer :: div, blockID, i, j, k, ngid
  logical :: allNegativeOneFound

  !-----------------------------------------------------------------------------
  ! build the tree information from the gid array and the number of blocks
  ! on each processor
  !-----------------------------------------------------------------------------

  div = xx*alnblocks

  allNegativeOneFound = .FALSE.
  do blockID = 1,localNumBlocks

     ! neighbor data
     ngid = 0
     do j = 1,nfaces
        ngid = ngid + 1

        if (gr_gid(ngid,blockID).gt.0) then

           if (gr_gid(ngid,blockID).le.div) then
              neigh(2,j,blockID) = int((gr_gid(ngid,blockID)-1)/alnblocks)

              if (neigh(2,j,blockID).gt.gr_meshNumProcs-1)  &
                   neigh(2,j,blockID) = gr_meshNumProcs - 1

              neigh(1,j,blockID) = gr_gid(ngid,blockID) -  &
                   (alnblocks*neigh(2,j,blockID))
           else
              neigh(2,j,blockID) = &
                   int((gr_gid(ngid,blockID)-1-div)/(alnblocks-1)) + xx

              if (neigh(2,j,blockID).gt.gr_meshNumProcs-1)  &
                   neigh(2,j,blockID) = gr_meshNumProcs - 1

              neigh(1,j,blockID) = gr_gid(ngid,blockID) - div - &
                   ((alnblocks-1)*(neigh(2,j,blockID)-xx))
           end if
        else
           neigh(1,j,blockID) = gr_gid(ngid,blockID)
           neigh(2,j,blockID) = gr_gid(ngid,blockID)
        end if
     end do

     ! parent data
     ngid = ngid + 1
     if (gr_gid(ngid,blockID).gt.0) then
        if (gr_gid(ngid,blockID).le.div) then
           parent(2,blockID) = int((gr_gid(ngid,blockID)-1)/alnblocks)
           if (parent(2,blockID).gt.gr_meshNumProcs-1)  &
                parent(2,blockID) = gr_meshNumProcs - 1
           parent(1,blockID) = gr_gid(ngid,blockID) -  &
                (alnblocks*parent(2,blockID))
        else
           parent(2,blockID) = &
                int((gr_gid(ngid,blockID)-1-div)/(alnblocks-1)) + xx
           if (parent(2,blockID).gt.gr_meshNumProcs-1)  &
                parent(2,blockID) = gr_meshNumProcs - 1
           parent(1,blockID) = gr_gid(ngid,blockID) - div - &
                ((alnblocks-1)*(parent(2,blockID)-xx))
        end if
     else
        parent(1,blockID) = gr_gid(ngid,blockID)
        parent(2,blockID) = gr_gid(ngid,blockID)
     end if

     ! children data
     do j = 1,nchild
        ngid = ngid + 1
        if (gr_gid(ngid,blockID).gt.0) then
           if (gr_gid(ngid,blockID).le.div) then
              child(2,j,blockID) = int((gr_gid(ngid,blockID)-1)/alnblocks)
              if (child(2,j,blockID).gt.gr_meshNumProcs-1)  &
                   child(2,j,blockID) = gr_meshNumProcs - 1
              child(1,j,blockID) = gr_gid(ngid,blockID) -  &
                   (alnblocks*child(2,j,blockID))
           else
              child(2,j,blockID) = &
                   int((gr_gid(ngid,blockID)-1-div)/(alnblocks-1)) + xx
              if (child(2,j,blockID).gt.gr_meshNumProcs-1)  &
                   child(2,j,blockID) = gr_meshNumProcs - 1
              child(1,j,blockID) = gr_gid(ngid,blockID) - div - &
                   ((alnblocks-1)*(child(2,j,blockID)-xx))
           end if
        else
           child(1,j,blockID) = gr_gid(ngid,blockID)
           child(2,j,blockID) = gr_gid(ngid,blockID)
        end if
     end do

     ! If ALL of the gr_gid information that we just read in and then
     ! used to initialize the neigh, parent, and child arrays for a block
     ! consists of '-1' values, then this tree metainformation cannot
     ! be part of the description of a valid Paramesh4 Grid
     ! configuration. [It can happen in a valid grid configuration that
     ! all of neigh for a block is -1, and it can happen that all of
     ! parent for a block is -1 (if the parent is a root block); but
     ! both together cannot validly happen for the same block. If a
     ! root blocks were to have no same-level neighbors according to
     ! neigh, i.e., neigh filled with negative values, then those
     ! values should be <= -20 indicating boundary conditions rather
     ! than -1.]
     ! We may be restarting from a checkpoint that was written by
     ! FLASH-with-Amrex.
     if (ALL(gr_gid(1:ngid,blockID) == -1)) allNegativeOneFound = .TRUE.
  end do

  gr_gidIsValid = .NOT. allNegativeOneFound

  gsurrblks_set = -1

#ifdef FLASH_GRID_PARAMESH4DEV_SURR_BLKS_OPTION
  ! surrounding neighbor data
  if (gr_is_gsurr_blks_initialized) then
     if (gr_globalMe == MASTER_PE) &
          print *, "read 'gsurr_blks' dataset - applying pm4dev optimization."

     allNegativeOneFound = .FALSE.
     do blockID = 1,localNumBlocks
        do k = 1, 1+(K3D*2)
           do j = 1, 1+(K2D*2)
              do i = 1, 1+(K1D*2)
                 !Store the node type first of all.
                 surr_blks(3,i,j,k,blockID) = &
                      gr_gsurr_blks(2,i,j,k,blockID)
                 if (gr_gsurr_blks(1,i,j,k,blockID).gt.0) then
                    if (gr_gsurr_blks(1,i,j,k,blockID).le.div) then
                       surr_blks(2,i,j,k,blockID) = &
                            int((gr_gsurr_blks(1,i,j,k,blockID)-1)/alnblocks)
                       if (surr_blks(2,i,j,k,blockID).gt.gr_meshNumProcs-1)  &
                            surr_blks(2,i,j,k,blockID) = gr_meshNumProcs - 1
                       surr_blks(1,i,j,k,blockID) = &
                            gr_gsurr_blks(1,i,j,k,blockID) - &
                            (alnblocks*surr_blks(2,i,j,k,blockID))
                    else
                       surr_blks(2,i,j,k,blockID) = &
                            int((gr_gsurr_blks(1,i,j,k,blockID)-1-div)/&
                            (alnblocks-1)) + xx
                       if (surr_blks(2,i,j,k,blockID).gt.gr_meshNumProcs-1)  &
                            surr_blks(2,i,j,k,blockID) = gr_meshNumProcs - 1
                       surr_blks(1,i,j,k,blockID) = &
                            gr_gsurr_blks(1,i,j,k,blockID) - div - &
                            ((alnblocks-1)*(surr_blks(2,i,j,k,blockID)-xx))
                    end if
                 else
                    surr_blks(1,i,j,k,blockID) = gr_gsurr_blks(1,i,j,k,blockID)
                    surr_blks(2,i,j,k,blockID) = gr_gsurr_blks(1,i,j,k,blockID)
                    surr_blks(3,i,j,k,blockID) = gr_gsurr_blks(2,i,j,k,blockID)
                 end if
              end do
           end do
        end do
        ! If ALL of the surr_blks information that we just filled in for a block
        ! consists of '-1' values, then surr_blks cannot be valid.
        ! We may be restarting from a checkpoint that was not written by
        ! FLASH-with-Paramesh4.
        if (ALL(surr_blks(1:3,1:3,1:1+(K2D*2),1:1+(K3D*2),blockID) == -1)) &
             allNegativeOneFound = .TRUE.
     end do

     if (.NOT. allNegativeOneFound) gsurrblks_set = +1

     if (use_flash_surr_blks_fill .AND. localNumBlocks > 0) then
        surr_blks_valid = .NOT. allNegativeOneFound
     else
        surr_blks_valid = use_flash_surr_blks_fill
     end if
  else
     surr_blks_valid = .false.
  end if
#endif

#ifdef DEBUG_IO
  print*,' @g',gr_globalMe,' gr_gidIsValid=',gr_gidIsValid,&
                           ' localNumBlocks=',localNumBlocks,&
                           ' gr_is_gsurr_blks_initialized=',gr_is_gsurr_blks_initialized,&
                           ' use_flash_surr_blks_fill=',use_flash_surr_blks_fill,&
                           ' surr_blks_valid=',surr_blks_valid,&
                           ' gsurrblks_set=',gsurrblks_set
#endif
end subroutine Grid_receiveInputData
