!!****if* source/numericalTools/MoL/MoLMemory/ml_memAddToVarsImpl
!!
!!  NAME
!!
!!      ml_memAddToVarsImpl
!!
!!  SYNOPSIS
!!
!!      use ml_memAddToVarsImpl, only:
!!
!!  DESCRIPTION
!!
!!      Perform a linear combination of source terms into the specified destination
!!      for all variables evolved by MoL.  The destination can either be the evolved
!!      variables in UNK or their corresponding locations in one of the MoL-specific
!!      scratch memory locations.  The source terms will only-ever be taken from
!!      MoL-specific scratch memory locations.  The linear combinations will take one
!!      of the following forms:
!!
!!          dst = dstFac*dst + val
!!
!!              or
!!
!!          dst = dstFac*dst + fac1*src1 + ... + facN*srcN
!!
!!      Valid source and destination locations include (defined in MoL.h):
!!          - MOL_EVOLVED : Evolved variables in UNK
!!          - MOL_INITIAL : Copy of the evolved variables at the start of a timestep
!!          - MOL_RHS     : The currently-being-calculated RHS terms
!!          - other       : Each integrator may specify some additional number of
!!                          of scratch-memory for intermediate stages/RHS terms
!!
!!  ARGUMENTS
!!
!!      TODO: change to doxygen
!!
!!  NOTES
!!
!!      For optimal memory use, there are specific implementations for up to twelve
!!      source terms - this will need to increase if any methods beyond 4th-order
!!      IMEX(-MRI) are added in the future.  For now, twelve source terms are
!!      sufficient for the 4rd-order IMEX-MRI-GARK multi-rate integrator
!!***
!!REORDER(4): dstPtr, src[12345678]Ptr
!!REORDER(4): src10Ptr, src11Ptr, src12Ptr
module ml_memAddToVarsImpl

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

   implicit none

   interface ml_memAddToVarsN
      procedure :: ml_memAddToVars0, ml_memAddToVars1, ml_memAddToVars2, &
         ml_memAddToVars3, ml_memAddToVars4, ml_memAddToVars5, &
         ml_memAddToVars6, ml_memAddToVars7, ml_memAddToVars8, &
         ml_memAddToVars9, ml_memAddToVars10, ml_memAddToVars11, ml_memAddToVars12
   end interface ml_memAddToVarsN

contains

   subroutine ml_memAddToVars0(dst, dstFac, val)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      real, intent(in) :: val

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
               do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                  dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) + val
               end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars0

   subroutine ml_memAddToVars1(dst, dstFac, src1, fac1)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      integer, intent(in) :: src1
      real, intent(in) :: fac1

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr
      real, dimension(:, :, :, :), pointer :: src1Ptr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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
         call MoL_getDataPtr(tileDesc, src1Ptr, src1)

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
               do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                  dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) &
                                             + fac1*src1Ptr(:, i, j, k)
               end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, src1Ptr, src1)
         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars1

   subroutine ml_memAddToVars2(dst, dstFac, src1, src2, fac1, fac2)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      integer, intent(in) :: src1, src2
      real, intent(in) :: fac1, fac2

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr
      real, dimension(:, :, :, :), pointer :: src1Ptr, src2Ptr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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
         call MoL_getDataPtr(tileDesc, src1Ptr, src1)
         call MoL_getDataPtr(tileDesc, src2Ptr, src2)

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
               do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                  dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) &
                                             + fac1*src1Ptr(:, i, j, k) &
                                             + fac2*src2Ptr(:, i, j, k)
               end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, src2Ptr, src2)
         call MoL_releaseDataPtr(tileDesc, src1Ptr, src1)
         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars2

   subroutine ml_memAddToVars3(dst, dstFac, src1, src2, src3, &
                               fac1, fac2, fac3)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      integer, intent(in) :: src1, src2, src3
      real, intent(in) :: fac1, fac2, fac3

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr
      real, dimension(:, :, :, :), pointer :: src1Ptr, src2Ptr, src3Ptr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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
         call MoL_getDataPtr(tileDesc, src1Ptr, src1)
         call MoL_getDataPtr(tileDesc, src2Ptr, src2)
         call MoL_getDataPtr(tileDesc, src3Ptr, src3)

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
               do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                  dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) &
                                             + fac1*src1Ptr(:, i, j, k) &
                                             + fac2*src2Ptr(:, i, j, k) &
                                             + fac3*src3Ptr(:, i, j, k)
               end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, src3Ptr, src3)
         call MoL_releaseDataPtr(tileDesc, src2Ptr, src2)
         call MoL_releaseDataPtr(tileDesc, src1Ptr, src1)
         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars3

   subroutine ml_memAddToVars4(dst, dstFac, src1, src2, src3, src4, &
                               fac1, fac2, fac3, fac4)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      integer, intent(in) :: src1, src2, src3, src4
      real, intent(in) :: fac1, fac2, fac3, fac4

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr
      real, dimension(:, :, :, :), pointer :: src1Ptr, src2Ptr, src3Ptr, src4Ptr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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
         call MoL_getDataPtr(tileDesc, src1Ptr, src1)
         call MoL_getDataPtr(tileDesc, src2Ptr, src2)
         call MoL_getDataPtr(tileDesc, src3Ptr, src3)
         call MoL_getDataPtr(tileDesc, src4Ptr, src4)

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
               do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                  dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) &
                                             + fac1*src1Ptr(:, i, j, k) &
                                             + fac2*src2Ptr(:, i, j, k) &
                                             + fac3*src3Ptr(:, i, j, k) &
                                             + fac4*src4Ptr(:, i, j, k)
               end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, src4Ptr, src4)
         call MoL_releaseDataPtr(tileDesc, src3Ptr, src3)
         call MoL_releaseDataPtr(tileDesc, src2Ptr, src2)
         call MoL_releaseDataPtr(tileDesc, src1Ptr, src1)
         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars4

   subroutine ml_memAddToVars5(dst, dstFac, src1, src2, src3, src4, src5, &
                               fac1, fac2, fac3, fac4, fac5)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      integer, intent(in) :: src1, src2, src3, src4, src5
      real, intent(in) :: fac1, fac2, fac3, fac4, fac5

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr
      real, dimension(:, :, :, :), pointer :: src1Ptr, src2Ptr, src3Ptr, src4Ptr, src5Ptr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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
         call MoL_getDataPtr(tileDesc, src1Ptr, src1)
         call MoL_getDataPtr(tileDesc, src2Ptr, src2)
         call MoL_getDataPtr(tileDesc, src3Ptr, src3)
         call MoL_getDataPtr(tileDesc, src4Ptr, src4)
         call MoL_getDataPtr(tileDesc, src5Ptr, src5)

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
               do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                  dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) &
                                             + fac1*src1Ptr(:, i, j, k) &
                                             + fac2*src2Ptr(:, i, j, k) &
                                             + fac3*src3Ptr(:, i, j, k) &
                                             + fac4*src4Ptr(:, i, j, k) &
                                             + fac5*src5Ptr(:, i, j, k)
               end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, src5Ptr, src5)
         call MoL_releaseDataPtr(tileDesc, src4Ptr, src4)
         call MoL_releaseDataPtr(tileDesc, src3Ptr, src3)
         call MoL_releaseDataPtr(tileDesc, src2Ptr, src2)
         call MoL_releaseDataPtr(tileDesc, src1Ptr, src1)
         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars5

   subroutine ml_memAddToVars6(dst, dstFac, src1, src2, src3, src4, src5, src6, &
                               fac1, fac2, fac3, fac4, fac5, fac6)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      integer, intent(in) :: src1, src2, src3, src4, src5, src6
      real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr
      real, dimension(:, :, :, :), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                              src4Ptr, src5Ptr, src6Ptr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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
         call MoL_getDataPtr(tileDesc, src1Ptr, src1)
         call MoL_getDataPtr(tileDesc, src2Ptr, src2)
         call MoL_getDataPtr(tileDesc, src3Ptr, src3)
         call MoL_getDataPtr(tileDesc, src4Ptr, src4)
         call MoL_getDataPtr(tileDesc, src5Ptr, src5)
         call MoL_getDataPtr(tileDesc, src6Ptr, src6)

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
               do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                  dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) &
                                             + fac1*src1Ptr(:, i, j, k) &
                                             + fac2*src2Ptr(:, i, j, k) &
                                             + fac3*src3Ptr(:, i, j, k) &
                                             + fac4*src4Ptr(:, i, j, k) &
                                             + fac5*src5Ptr(:, i, j, k) &
                                             + fac6*src6Ptr(:, i, j, k)
               end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, src6Ptr, src6)
         call MoL_releaseDataPtr(tileDesc, src5Ptr, src5)
         call MoL_releaseDataPtr(tileDesc, src4Ptr, src4)
         call MoL_releaseDataPtr(tileDesc, src3Ptr, src3)
         call MoL_releaseDataPtr(tileDesc, src2Ptr, src2)
         call MoL_releaseDataPtr(tileDesc, src1Ptr, src1)
         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars6

   subroutine ml_memAddToVars7(dst, dstFac, &
                               src1, src2, src3, src4, src5, src6, src7, &
                               fac1, fac2, fac3, fac4, fac5, fac6, fac7)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      integer, intent(in) :: src1, src2, src3, src4, src5, src6, src7
      real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr
      real, dimension(:, :, :, :), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                              src4Ptr, src5Ptr, src6Ptr, &
                                              src7Ptr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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
         call MoL_getDataPtr(tileDesc, src1Ptr, src1)
         call MoL_getDataPtr(tileDesc, src2Ptr, src2)
         call MoL_getDataPtr(tileDesc, src3Ptr, src3)
         call MoL_getDataPtr(tileDesc, src4Ptr, src4)
         call MoL_getDataPtr(tileDesc, src5Ptr, src5)
         call MoL_getDataPtr(tileDesc, src6Ptr, src6)
         call MoL_getDataPtr(tileDesc, src7Ptr, src7)

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
            do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
               dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) &
                                          + fac1*src1Ptr(:, i, j, k) &
                                          + fac2*src2Ptr(:, i, j, k) &
                                          + fac3*src3Ptr(:, i, j, k) &
                                          + fac4*src4Ptr(:, i, j, k) &
                                          + fac5*src5Ptr(:, i, j, k) &
                                          + fac6*src6Ptr(:, i, j, k) &
                                          + fac7*src7Ptr(:, i, j, k)
            end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, src7Ptr, src7)
         call MoL_releaseDataPtr(tileDesc, src6Ptr, src6)
         call MoL_releaseDataPtr(tileDesc, src5Ptr, src5)
         call MoL_releaseDataPtr(tileDesc, src4Ptr, src4)
         call MoL_releaseDataPtr(tileDesc, src3Ptr, src3)
         call MoL_releaseDataPtr(tileDesc, src2Ptr, src2)
         call MoL_releaseDataPtr(tileDesc, src1Ptr, src1)
         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars7

   subroutine ml_memAddToVars8(dst, dstFac, &
                               src1, src2, src3, src4, src5, src6, src7, src8, &
                               fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      integer, intent(in) :: src1, src2, src3, src4, src5, src6, src7, src8
      real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr
      real, dimension(:, :, :, :), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                              src4Ptr, src5Ptr, src6Ptr, &
                                              src7Ptr, src8Ptr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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
         call MoL_getDataPtr(tileDesc, src1Ptr, src1)
         call MoL_getDataPtr(tileDesc, src2Ptr, src2)
         call MoL_getDataPtr(tileDesc, src3Ptr, src3)
         call MoL_getDataPtr(tileDesc, src4Ptr, src4)
         call MoL_getDataPtr(tileDesc, src5Ptr, src5)
         call MoL_getDataPtr(tileDesc, src6Ptr, src6)
         call MoL_getDataPtr(tileDesc, src7Ptr, src7)
         call MoL_getDataPtr(tileDesc, src8Ptr, src8)

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
               do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                  dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) &
                                             + fac1*src1Ptr(:, i, j, k) &
                                             + fac2*src2Ptr(:, i, j, k) &
                                             + fac3*src3Ptr(:, i, j, k) &
                                             + fac4*src4Ptr(:, i, j, k) &
                                             + fac5*src5Ptr(:, i, j, k) &
                                             + fac6*src6Ptr(:, i, j, k) &
                                             + fac7*src7Ptr(:, i, j, k) &
                                             + fac8*src8Ptr(:, i, j, k)
               end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, src8Ptr, src8)
         call MoL_releaseDataPtr(tileDesc, src7Ptr, src7)
         call MoL_releaseDataPtr(tileDesc, src6Ptr, src6)
         call MoL_releaseDataPtr(tileDesc, src5Ptr, src5)
         call MoL_releaseDataPtr(tileDesc, src4Ptr, src4)
         call MoL_releaseDataPtr(tileDesc, src3Ptr, src3)
         call MoL_releaseDataPtr(tileDesc, src2Ptr, src2)
         call MoL_releaseDataPtr(tileDesc, src1Ptr, src1)
         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars8

   subroutine ml_memAddToVars9(dst, dstFac, &
                               src1, src2, src3, src4, src5, src6, src7, src8, src9, &
                               fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      integer, intent(in) :: src1, src2, src3, src4, src5, src6, src7, src8, src9
      real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr
      real, dimension(:, :, :, :), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                              src4Ptr, src5Ptr, src6Ptr, &
                                              src7Ptr, src8Ptr, src9Ptr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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
         call MoL_getDataPtr(tileDesc, src1Ptr, src1)
         call MoL_getDataPtr(tileDesc, src2Ptr, src2)
         call MoL_getDataPtr(tileDesc, src3Ptr, src3)
         call MoL_getDataPtr(tileDesc, src4Ptr, src4)
         call MoL_getDataPtr(tileDesc, src5Ptr, src5)
         call MoL_getDataPtr(tileDesc, src6Ptr, src6)
         call MoL_getDataPtr(tileDesc, src7Ptr, src7)
         call MoL_getDataPtr(tileDesc, src8Ptr, src8)
         call MoL_getDataPtr(tileDesc, src9Ptr, src9)

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
               do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                  dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) &
                                             + fac1*src1Ptr(:, i, j, k) &
                                             + fac2*src2Ptr(:, i, j, k) &
                                             + fac3*src3Ptr(:, i, j, k) &
                                             + fac4*src4Ptr(:, i, j, k) &
                                             + fac5*src5Ptr(:, i, j, k) &
                                             + fac6*src6Ptr(:, i, j, k) &
                                             + fac7*src7Ptr(:, i, j, k) &
                                             + fac8*src8Ptr(:, i, j, k) &
                                             + fac9*src9Ptr(:, i, j, k)
               end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, src9Ptr, src9)
         call MoL_releaseDataPtr(tileDesc, src8Ptr, src8)
         call MoL_releaseDataPtr(tileDesc, src7Ptr, src7)
         call MoL_releaseDataPtr(tileDesc, src6Ptr, src6)
         call MoL_releaseDataPtr(tileDesc, src5Ptr, src5)
         call MoL_releaseDataPtr(tileDesc, src4Ptr, src4)
         call MoL_releaseDataPtr(tileDesc, src3Ptr, src3)
         call MoL_releaseDataPtr(tileDesc, src2Ptr, src2)
         call MoL_releaseDataPtr(tileDesc, src1Ptr, src1)
         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars9

   subroutine ml_memAddToVars10(dst, dstFac, &
                                src1, src2, src3, src4, src5, src6, src7, src8, &
                                src9, src10, &
                                fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, &
                                fac9, fac10)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      integer, intent(in) :: src1, src2, src3, src4, src5, src6, src7, src8, src9, src10
      real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9, fac10

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr
      real, dimension(:, :, :, :), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                              src4Ptr, src5Ptr, src6Ptr, &
                                              src7Ptr, src8Ptr, src9Ptr, &
                                              src10Ptr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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
         call MoL_getDataPtr(tileDesc, src1Ptr, src1)
         call MoL_getDataPtr(tileDesc, src2Ptr, src2)
         call MoL_getDataPtr(tileDesc, src3Ptr, src3)
         call MoL_getDataPtr(tileDesc, src4Ptr, src4)
         call MoL_getDataPtr(tileDesc, src5Ptr, src5)
         call MoL_getDataPtr(tileDesc, src6Ptr, src6)
         call MoL_getDataPtr(tileDesc, src7Ptr, src7)
         call MoL_getDataPtr(tileDesc, src8Ptr, src8)
         call MoL_getDataPtr(tileDesc, src9Ptr, src9)
         call MoL_getDataPtr(tileDesc, src10Ptr, src10)

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
               do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                  dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) &
                                             + fac1*src1Ptr(:, i, j, k) &
                                             + fac2*src2Ptr(:, i, j, k) &
                                             + fac3*src3Ptr(:, i, j, k) &
                                             + fac4*src4Ptr(:, i, j, k) &
                                             + fac5*src5Ptr(:, i, j, k) &
                                             + fac6*src6Ptr(:, i, j, k) &
                                             + fac7*src7Ptr(:, i, j, k) &
                                             + fac8*src8Ptr(:, i, j, k) &
                                             + fac9*src9Ptr(:, i, j, k) &
                                             + fac10*src10Ptr(:, i, j, k)
               end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, src10Ptr, src10)
         call MoL_releaseDataPtr(tileDesc, src9Ptr, src9)
         call MoL_releaseDataPtr(tileDesc, src8Ptr, src8)
         call MoL_releaseDataPtr(tileDesc, src7Ptr, src7)
         call MoL_releaseDataPtr(tileDesc, src6Ptr, src6)
         call MoL_releaseDataPtr(tileDesc, src5Ptr, src5)
         call MoL_releaseDataPtr(tileDesc, src4Ptr, src4)
         call MoL_releaseDataPtr(tileDesc, src3Ptr, src3)
         call MoL_releaseDataPtr(tileDesc, src2Ptr, src2)
         call MoL_releaseDataPtr(tileDesc, src1Ptr, src1)
         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars10

   subroutine ml_memAddToVars11(dst, dstFac, &
                                src1, src2, src3, src4, src5, src6, src7, src8, &
                                src9, src10, src11, &
                                fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, &
                                fac9, fac10, fac11)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                             src7, src8, src9, src10, src11
      real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                          fac7, fac8, fac9, fac10, fac11

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr
      real, dimension(:, :, :, :), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                              src4Ptr, src5Ptr, src6Ptr, &
                                              src7Ptr, src8Ptr, src9Ptr, &
                                              src10Ptr, src11Ptr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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
         call MoL_getDataPtr(tileDesc, src1Ptr, src1)
         call MoL_getDataPtr(tileDesc, src2Ptr, src2)
         call MoL_getDataPtr(tileDesc, src3Ptr, src3)
         call MoL_getDataPtr(tileDesc, src4Ptr, src4)
         call MoL_getDataPtr(tileDesc, src5Ptr, src5)
         call MoL_getDataPtr(tileDesc, src6Ptr, src6)
         call MoL_getDataPtr(tileDesc, src7Ptr, src7)
         call MoL_getDataPtr(tileDesc, src8Ptr, src8)
         call MoL_getDataPtr(tileDesc, src9Ptr, src9)
         call MoL_getDataPtr(tileDesc, src10Ptr, src10)
         call MoL_getDataPtr(tileDesc, src11Ptr, src11)

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
               do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                  dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) &
                                             + fac1*src1Ptr(:, i, j, k) &
                                             + fac2*src2Ptr(:, i, j, k) &
                                             + fac3*src3Ptr(:, i, j, k) &
                                             + fac4*src4Ptr(:, i, j, k) &
                                             + fac5*src5Ptr(:, i, j, k) &
                                             + fac6*src6Ptr(:, i, j, k) &
                                             + fac7*src7Ptr(:, i, j, k) &
                                             + fac8*src8Ptr(:, i, j, k) &
                                             + fac9*src9Ptr(:, i, j, k) &
                                             + fac10*src10Ptr(:, i, j, k) &
                                             + fac11*src11Ptr(:, i, j, k)
               end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, src11Ptr, src11)
         call MoL_releaseDataPtr(tileDesc, src10Ptr, src10)
         call MoL_releaseDataPtr(tileDesc, src9Ptr, src9)
         call MoL_releaseDataPtr(tileDesc, src8Ptr, src8)
         call MoL_releaseDataPtr(tileDesc, src7Ptr, src7)
         call MoL_releaseDataPtr(tileDesc, src6Ptr, src6)
         call MoL_releaseDataPtr(tileDesc, src5Ptr, src5)
         call MoL_releaseDataPtr(tileDesc, src4Ptr, src4)
         call MoL_releaseDataPtr(tileDesc, src3Ptr, src3)
         call MoL_releaseDataPtr(tileDesc, src2Ptr, src2)
         call MoL_releaseDataPtr(tileDesc, src1Ptr, src1)
         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars11

   subroutine ml_memAddToVars12(dst, dstFac, &
                                src1, src2, src3, src4, src5, src6, src7, src8, &
                                src9, src10, src11, src12, &
                                fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, &
                                fac9, fac10, fac11, fac12)
      use ml_variables, only: ml_nvars, ml_unk_mask, ml_scratch_mask
      use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

      use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
      use Grid_tile, only: Grid_tile_t
      use Grid_iterator, only: Grid_iterator_t

      implicit none

      integer, intent(in) :: dst
      real, intent(in) :: dstFac
      integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                             src7, src8, src9, src10, src11, src12
      real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                          fac7, fac8, fac9, fac10, fac11, fac12

      integer, dimension(ml_nvars) :: dstVars
      real, dimension(:, :, :, :), pointer :: dstPtr
      real, dimension(:, :, :, :), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                              src4Ptr, src5Ptr, src6Ptr, &
                                              src7Ptr, src8Ptr, src9Ptr, &
                                              src10Ptr, src11Ptr, src12Ptr

      type(Grid_iterator_t) :: itor
      type(Grid_tile_t) :: tileDesc

      integer :: i, j, k

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
         call MoL_getDataPtr(tileDesc, src1Ptr, src1)
         call MoL_getDataPtr(tileDesc, src2Ptr, src2)
         call MoL_getDataPtr(tileDesc, src3Ptr, src3)
         call MoL_getDataPtr(tileDesc, src4Ptr, src4)
         call MoL_getDataPtr(tileDesc, src5Ptr, src5)
         call MoL_getDataPtr(tileDesc, src6Ptr, src6)
         call MoL_getDataPtr(tileDesc, src7Ptr, src7)
         call MoL_getDataPtr(tileDesc, src8Ptr, src8)
         call MoL_getDataPtr(tileDesc, src9Ptr, src9)
         call MoL_getDataPtr(tileDesc, src10Ptr, src10)
         call MoL_getDataPtr(tileDesc, src11Ptr, src11)
         call MoL_getDataPtr(tileDesc, src12Ptr, src12)

         do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
            do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
               do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                  dstPtr(dstVars, i, j, k) = dstFac*dstPtr(dstVars, i, j, k) &
                                             + fac1*src1Ptr(:, i, j, k) &
                                             + fac2*src2Ptr(:, i, j, k) &
                                             + fac3*src3Ptr(:, i, j, k) &
                                             + fac4*src4Ptr(:, i, j, k) &
                                             + fac5*src5Ptr(:, i, j, k) &
                                             + fac6*src6Ptr(:, i, j, k) &
                                             + fac7*src7Ptr(:, i, j, k) &
                                             + fac8*src8Ptr(:, i, j, k) &
                                             + fac9*src9Ptr(:, i, j, k) &
                                             + fac10*src10Ptr(:, i, j, k) &
                                             + fac11*src11Ptr(:, i, j, k) &
                                             + fac12*src12Ptr(:, i, j, k)
               end do ! i
            end do ! j
         end do ! k

         call MoL_releaseDataPtr(tileDesc, src12Ptr, src12)
         call MoL_releaseDataPtr(tileDesc, src11Ptr, src11)
         call MoL_releaseDataPtr(tileDesc, src10Ptr, src10)
         call MoL_releaseDataPtr(tileDesc, src9Ptr, src9)
         call MoL_releaseDataPtr(tileDesc, src8Ptr, src8)
         call MoL_releaseDataPtr(tileDesc, src7Ptr, src7)
         call MoL_releaseDataPtr(tileDesc, src6Ptr, src6)
         call MoL_releaseDataPtr(tileDesc, src5Ptr, src5)
         call MoL_releaseDataPtr(tileDesc, src4Ptr, src4)
         call MoL_releaseDataPtr(tileDesc, src3Ptr, src3)
         call MoL_releaseDataPtr(tileDesc, src2Ptr, src2)
         call MoL_releaseDataPtr(tileDesc, src1Ptr, src1)
         call MoL_releaseDataPtr(tileDesc, dstPtr, dst)

         call itor%next()
      end do TileLoop

      call Grid_releaseTileIterator(itor)
   end subroutine ml_memAddToVars12

end module ml_memAddToVarsImpl
