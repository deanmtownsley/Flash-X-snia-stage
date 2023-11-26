!!****if* source/Grid/GridMain/AMR/Amrex/gr_leafBlockInfo
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  gr_leafBlockInfo
!!
!! SYNOPSIS
!!
!!  use gr_leafBlockInfo
!!
!!  call gr_leafBlockInfoUpdate()
!!  leafID = gr_getLeafBlockNo(Grid_tile_t(IN) :: blockDesc)
!!
!! DESCRIPTION
!!
!!  This module provides a mapping from
!!       pairs (level, gid), representing
!!       refinement level and AMReX grid_index within the Grid multifab for
!!       that level, and together identifying a block globally and
!!       completely,
!!  to
!!       a local leaf ID, representing uniquely each relevant LEAF block that
!!       belongs to the executing MPI task.
!!
!!  Leaf IDs are small positive integers that are assigned to local blocks
!!  here when gr_leafBlockInfoUpdate() is called.
!!
!!  Leaf IDs are useful for indexing into arrays (not managed by the
!!  underlying mesh implementation like AMReX) that hold various kinds of
!!  per-LEAF-block information.
!!
!! ARGUMENTS
!!
!!   blockDesc - block descriptor that holds level and grid_index of the
!!               block whose ID is sought.
!!
!! NOTES
!!
!!   gr_leafBlockInfoUpdate() should be called during program initialization
!!   (after the domain is initialized and the mesh has been organized and
!!   distributed), and then once each time after the Grid refinement is updated.
!!
!! SEE ALSO
!!
!!  Grid_putFluxData_block
!!  Grid_getFluxCorrData_block
!!  gr_auxFluxData
!!
!! HISTORY
!!  2023-10-29 K. Weide  created
!!
!!
!!***

#include "constants.h"
#include "Simulation.h"
module gr_leafBlockInfo
   use Grid_tile,           ONLY : Grid_tile_t
   implicit none

   integer,save,dimension(MAXBLOCKS) :: leafNoToGlobalIdx
   integer,save,dimension(MAXBLOCKS) :: leafNoToLevel
   integer,save :: numLeafBlocks ! Number of leaf nodes on this proc

contains
   subroutine gr_leafBlockInfoUpdate()
      use Driver_interface,    ONLY : Driver_abort
      use Grid_interface,      ONLY : Grid_getTileIterator
      use Grid_iterator,       ONLY : Grid_iterator_t
      implicit none
      type(Grid_iterator_t) :: itor
      type(Grid_tile_t)     :: tileDesc
      call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
      numLeafBlocks = 0
      do while(itor%isValid())
         call itor%currentTile(tileDesc)
         if (numLeafBlocks > 0) then
            if (tileDesc%level < leafNoToLevel(numLeafBlocks)) then
               call Driver_abort("gr_leafBlockInfoUpdate: levels not ascending!")
            else if (tileDesc%level == leafNoToLevel(numLeafBlocks) .AND. &
                     tileDesc%grid_index .LE. leafNoToGlobalIdx(numLeafBlocks)) then
               call Driver_abort("gr_leafBlockInfoUpdate: grid_index not ascending!")
            end if
         end if
         numLeafBlocks = numLeafBlocks + 1
         leafNoToGlobalIdx(numLeafBlocks) = tileDesc%grid_index
         leafNoToLevel(numLeafBlocks)     = tileDesc%level
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
   end subroutine gr_leafBlockInfoUpdate

   integer function gr_getLeafBlockNo(blockDesc)
      implicit none
      type(Grid_tile_t),intent(IN) :: blockDesc
      integer :: i
      gr_getLeafBlockNo = 0
      do i=1, numLeafBlocks
         if (leafNoToLevel(i) == blockDesc%level .AND. &
             leafNoToGlobalIdx(i) == blockDesc%grid_index) then
            gr_getLeafBlockNo = i
            return
         end if
      end do
   end function gr_getLeafBlockNo
end module gr_leafBlockInfo
! Local Variables:
! f90-program-indent: 3
! f90-associate-indent: 3
! indent-tabs-mode: nil
! End:
