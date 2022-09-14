!!****if* source/numericalTools/MoL/MoLMemory/ml_memAddToVars
!!
!!  NAME
!!
!!      ml_memAddToVars
!!
!!  SYNOPSIS
!!
!!      call ml_memAddToVars(integer, intent(in) :: dst
!!                           real,    intent(in) :: dstFac
!!                           integer, intent(in) :: nsrcs
!!                           integer, intent(in) :: srcs(nsrcs)
!!                           real,    intent(in) :: facs(nsrcs))
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
!!          dst = dstFac*dst + fac(1_*src(1_ + ... + fac(nsrcs)*src(nsrcs)
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
!!      dst    : Index of the destination location to store the linear combination
!!      dstFac : Scaling factor for the destination - set this to zero to overwrite
!!               the existing value
!!      val    : Scalar value to set all variables to in dst
!!      nsrcs  : Number of source terms for the general N-src implementation
!!      srcs   : Array of source index locations in MoL scratch memory
!!      facs   : Array of scaling factors for each source term
!!
!!  NOTES
!!
!!      For optimal memory use, there are specific implementations for up to twelve
!!      source terms - this will need to increase if any methods beyond 4th-order
!!      IMEX(-MRI) are added in the future.  For now, twelve source terms are
!!      sufficient for the 4rd-order IMEX-MRI-GARK multi-rate integrator
!!***
!!REORDER(4): dstPtr, srcPtr
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
