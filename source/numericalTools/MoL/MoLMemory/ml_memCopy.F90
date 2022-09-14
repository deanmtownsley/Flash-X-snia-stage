!!****if* source/numericalTools/MoL/MoLMemory/ml_memCopy
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
!!  NAME
!!
!!      ml_memCopy
!!
!!  SYNOPSIS
!!
!!      call ml_memCopy(integer, intent(in) :: dst
!!                      integer, intent(in) :: src)
!!
!!  DESCRIPTION
!!
!!      Copy the variables evolved by MoL from provided source to destination
!!      memory locations.  These memory locations either be UNK or MoL-specific
!!      scratch memory.
!!
!!      Valid locations include (defined in MoL.h):
!!          - MOL_EVOLVED : Evolved variables in UNK
!!          - MOL_INITIAL : Copy of the evolved variables at the start of a timestep
!!          - MOL_RHS     : The currently-being-calculated RHS terms
!!          - other       : Each integrator may specify some additional number of
!!                          of scratch-memory for intermediate stages/RHS terms
!!
!!  ARGUMENTS
!!
!!      dst : Index of the destination location to copy to
!!      src : Index of the source location to copy from
!!
!!***
!!REORDER(4): dstPtr, srcPtr
subroutine ml_memCopy(dst, src)
   use ml_variables, only: ml_scratch_mask, ml_unk_mask, ml_nvars
   use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
   use Grid_tile, only: Grid_tile_t
   use Grid_iterator, only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

   implicit none

   integer, intent(in) :: dst, src

   integer, dimension(ml_nvars) :: dstVars, srcVars
   real, dimension(:, :, :, :), pointer :: dstPtr, srcPtr

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   integer :: i, j, k

   ! Bail if we requested something stupid
   if (dst .eq. src) return

   ! Select the correct variable masks
   if (dst .eq. MOL_EVOLVED) then
      dstVars = ml_unk_mask
   else
      dstVars = ml_scratch_mask
   end if

   if (src .eq. MOL_EVOLVED) then
      srcVars = ml_unk_mask
   else
      srcVars = ml_scratch_mask
   end if

   call Grid_getTileIterator(itor, LEAF, tiling=.true.)

   TileLoop: do
      if (.not. itor%isValid()) exit TileLoop

      call itor%currentTile(tileDesc)

      call MoL_getDataPtr(tileDesc, dstPtr, dst)
      call MoL_getDataPtr(tileDesc, srcPtr, src)

      do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
         do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
            do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
               dstPtr(dstVars, i, j, k) = srcPtr(srcVars, i, j, k)
            end do ! i
         end do ! j
      end do ! k

      call MoL_releaseDataPtr(tileDesc, srcPtr, src)
      call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

      call itor%next()
   end do TileLoop

   call Grid_releaseTileIterator(itor)
end subroutine ml_memCopy
