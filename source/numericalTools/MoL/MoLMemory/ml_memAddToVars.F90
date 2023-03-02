!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief ml_memAddToVars implementation

!!REORDER(4): dstPtr, srcPtr

!> @ingroup MoLMemory
!!
!! @brief Implements ml_memAddToVars
!!
!! @stubref{ml_memAddToVars}
subroutine ml_memAddToVars(dst, dstFac, nsrcs, srcs, facs)
   use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
   use ml_memAddToVarsImpl, only: ml_memAddToVarsN
   use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
   use Grid_tile, only: Grid_tile_t
   use Grid_iterator, only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: nsrcs
   integer, intent(in) :: srcs(nsrcs)
   real, intent(in) :: facs(nsrcs)

   integer, dimension(ml_nvars) :: dstVars
   real, dimension(:, :, :, :), pointer :: dstPtr
   real, dimension(:, :, :, :), pointer :: srcPtr

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   integer :: i, j, k, n

   if (dst .eq. MOL_EVOLVED) then
      dstVars = ml_unk_mask
   else
      dstVars = ml_scratch_mask
   end if

   ! For now, just abort if nothing was passed in
   if (nsrcs .lt. 1) return

   select case (nsrcs)
   case (1)
      call ml_memAddToVarsN(dst, dstFac, &
                            srcs(1), &
                            facs(1))
   case (2)
      call ml_memAddToVarsN(dst, dstFac, &
                            srcs(1), srcs(2), &
                            facs(1), facs(2))
   case (3)
      call ml_memAddToVarsN(dst, dstFac, &
                            srcs(1), srcs(2), srcs(3), &
                            facs(1), facs(2), facs(3))

   case (4)
      call ml_memAddToVarsN(dst, dstFac, &
                            srcs(1), srcs(2), srcs(3), srcs(4), &
                            facs(1), facs(2), facs(3), facs(4))

   case (5)
      call ml_memAddToVarsN(dst, dstFac, &
                            srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                            facs(1), facs(2), facs(3), facs(4), facs(5))

   case (6)
      call ml_memAddToVarsN(dst, dstFac, &
                            srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), srcs(6), &
                            facs(1), facs(2), facs(3), facs(4), facs(5), facs(6))

   case (7)
      call ml_memAddToVarsN(dst, dstFac, &
                            srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                            srcs(6), srcs(7), &
                            facs(1), facs(2), facs(3), facs(4), facs(5), &
                            facs(6), facs(7))

   case (8)
      call ml_memAddToVarsN(dst, dstFac, &
                            srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                            srcs(6), srcs(7), srcs(8), &
                            facs(1), facs(2), facs(3), facs(4), facs(5), &
                            facs(6), facs(7), facs(8))

   case (9)
      call ml_memAddToVarsN(dst, dstFac, &
                            srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                            srcs(6), srcs(7), srcs(8), srcs(9), &
                            facs(1), facs(2), facs(3), facs(4), facs(5), &
                            facs(6), facs(7), facs(8), facs(9))

   case (10)
      call ml_memAddToVarsN(dst, dstFac, &
                            srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                            srcs(6), srcs(7), srcs(8), srcs(9), srcs(10), &
                            facs(1), facs(2), facs(3), facs(4), facs(5), &
                            facs(6), facs(7), facs(8), facs(9), facs(10))

   case (11)
      call ml_memAddToVarsN(dst, dstFac, &
                            srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                            srcs(6), srcs(7), srcs(8), srcs(9), srcs(10), &
                            srcs(11), &
                            facs(1), facs(2), facs(3), facs(4), facs(5), &
                            facs(6), facs(7), facs(8), facs(9), facs(10), &
                            facs(11))

   case (12)
      call ml_memAddToVarsN(dst, dstFac, &
                            srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                            srcs(6), srcs(7), srcs(8), srcs(9), srcs(10), &
                            srcs(11), srcs(12), &
                            facs(1), facs(2), facs(3), facs(4), facs(5), &
                            facs(6), facs(7), facs(8), facs(9), facs(10), &
                            facs(11), facs(12))

   case default
      ! Slow, but for now it will work
      if (dst .eq. MOL_EVOLVED) then
         dstVars = ml_unk_mask
      else
         dstVars = ml_scratch_mask
      end if

      call Grid_getTileIterator(itor, LEAF, tiling=.true.)

      TileLoop: do
         if (.not. itor%isValid()) exit TileLoop

         call itor%currentTile(tileDesc)

         call MoL_getDataPtr(tileDesc, dstPtr, dst)

         dstPtr(dstVars, :, :, :) = dstFac*dstPtr(dstVars, :, :, :)

         do n = 1, nsrcs
            call MoL_getDataPtr(tileDesc, srcPtr, srcs(n))

            do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
               do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
                  do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                     dstPtr(dstVars, i, j, k) = dstPtr(dstVars, i, j, k) &
                                                + facs(n)*srcPtr(:, i, j, k)
                  end do ! i
               end do ! j
            end do ! k

            call MoL_releaseDataPtr(tileDesc, srcPtr, srcs(n))
         end do ! n

         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end select
end subroutine ml_memAddToVars
