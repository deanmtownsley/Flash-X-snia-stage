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
!!      call ml_memAddToVars(integer, intent(in) :: dst
!!                           real,    intent(in) :: dstFac
!!                           real,    intent(in) :: val)
!!
!!      call ml_memAddToVars(integer, intent(in) :: dst
!!                           real,    intent(in) :: dstFac
!!                           integer, intent(in) :: src1, src2, ...
!!                           real,    intent(in) :: fac1, fac2, ...)
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
!!      src1
!!       |     : Index of each source term location to use
!!      srcN
!!      fac1
!!       |     : Scaling factor for each source term location to use
!!      facN
!!
!!  NOTES
!!
!!      For optimal memory use, there are specific implementations for up to twelve
!!      source terms - this will need to increase if any methods beyond 4th-order
!!      IMEX(-MRI) are added in the future.  For now, twelve source terms are
!!      sufficient for the 4rd-order IMEX-MRI-GARK multi-rate integrator
!!***
!!REORDER(4): dstPtr, srcPtr, src[12345678]Ptr
!!REORDER(4): src10Ptr, src11Ptr, src12Ptr
subroutine ml_memAddToVarsN(dst, dstFac, nsrcs, srcs, facs)
    use MoL_variables,  only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: nsrcs
    integer, intent(in) :: srcs(nsrcs)
    real,    intent(in) :: facs(nsrcs)

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: srcPtr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k, n

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    ! For now, just abort if nothing was passed in
    if (nsrcs .lt. 1) return

    select case(nsrcs)
    case (1)
        call ml_memAddToVars1(dst, dstFac, &
                               srcs(1),     &
                               facs(1))
    case (2)
        call ml_memAddToVars2(dst, dstFac,      &
                               srcs(1), srcs(2), &
                               facs(1), facs(2))
    case (3)
        call ml_memAddToVars3(dst, dstFac,               &
                               srcs(1), srcs(2), srcs(3), &
                               facs(1), facs(2), facs(3))

    case (4)
        call ml_memAddToVars4(dst, dstFac,                        &
                               srcs(1), srcs(2), srcs(3), srcs(4), &
                               facs(1), facs(2), facs(3), facs(4))

    case (5)
        call ml_memAddToVars5(dst, dstFac,                                 &
                               srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                               facs(1), facs(2), facs(3), facs(4), facs(5))

    case (6)
        call ml_memAddToVars6(dst, dstFac,                                          &
                               srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), srcs(6), &
                               facs(1), facs(2), facs(3), facs(4), facs(5), facs(6))

    case (7)
        call ml_memAddToVars7(dst, dstFac,                                 &
                               srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                               srcs(6), srcs(7),                            &
                               facs(1), facs(2), facs(3), facs(4), facs(5), &
                               facs(6), facs(7))

    case (8)
        call ml_memAddToVars8(dst, dstFac,                                 &
                               srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                               srcs(6), srcs(7), srcs(8),                   &
                               facs(1), facs(2), facs(3), facs(4), facs(5), &
                               facs(6), facs(7), facs(8))

    case (9)
        call ml_memAddToVars9(dst, dstFac,                                 &
                               srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                               srcs(6), srcs(7), srcs(8), srcs(9),          &
                               facs(1), facs(2), facs(3), facs(4), facs(5), &
                               facs(6), facs(7), facs(8), facs(9))

    case (10)
        call ml_memAddToVars10(dst, dstFac,                                 &
                                srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                                srcs(6), srcs(7), srcs(8), srcs(9), srcs(10),&
                                facs(1), facs(2), facs(3), facs(4), facs(5), &
                                facs(6), facs(7), facs(8), facs(9), facs(10))

    case (11)
        call ml_memAddToVars11(dst, dstFac,                                 &
                                srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                                srcs(6), srcs(7), srcs(8), srcs(9), srcs(10),&
                                srcs(11),                                    &
                                facs(1), facs(2), facs(3), facs(4), facs(5), &
                                facs(6), facs(7), facs(8), facs(9), facs(10),&
                                facs(11))

    case (12)
        call ml_memAddToVars12(dst, dstFac,                                 &
                                srcs(1), srcs(2), srcs(3), srcs(4), srcs(5), &
                                srcs(6), srcs(7), srcs(8), srcs(9), srcs(10),&
                                srcs(11), srcs(12),                          &
                                facs(1), facs(2), facs(3), facs(4), facs(5), &
                                facs(6), facs(7), facs(8), facs(9), facs(10),&
                                facs(11), facs(12))

    case default
        ! Slow, but for now it will work
        if (dst .eq. MOL_EVOLVED) then
            dstVars = MoL_unk_mask
        else
            dstVars = MoL_scratch_mask
        end if
    
        call Grid_getTileIterator(itor, LEAF)
    
        TileLoop: do
            if (.not. itor%isValid()) exit TileLoop
    
            call itor%currentTile(tileDesc)
    
            call ml_memGetDataPtr(tileDesc, dstPtr,  dst)

            dstPtr(dstVars,:,:,:) = dstFac*dstPtr(dstVars,:,:,:)
            
            do n = 1, nsrcs
                call ml_memGetDataPtr(tileDesc, srcPtr, srcs(n))

                do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
                    do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                        do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                            dstPtr(dstVars,i,j,k) =           dstPtr(dstVars,i,j,k) &
                                                    + facs(n)*srcPtr(      :,i,j,k)
                        end do ! i
                    end do ! j
                end do ! k

                call ml_memReleaseDataPtr(tileDesc, srcPtr, srcs(n))
            end do ! n
    
            call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)
    
            call itor%next()
        end do TileLoop
    
        call Grid_releaseTileIterator(itor)
    end select
end subroutine ml_memAddToVarsN

subroutine ml_memAddToVars0(dst, dstFac, val)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    real,    intent(in) :: val

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr, dst)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
            do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) = dstFac*dstPtr(dstVars,i,j,k) + val
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, dstPtr, dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars0

subroutine ml_memAddToVars1(dst, dstFac, src1, fac1)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: src1
    real,    intent(in) :: fac1

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: src1Ptr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr,  dst)
        call ml_memGetDataPtr(tileDesc, src1Ptr, src1)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
            do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) =  dstFac*dstPtr(dstVars,i,j,k) &
                                            + fac1*src1Ptr(      :,i,j,k)
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, src1Ptr, src1)
        call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars1

subroutine ml_memAddToVars2(dst, dstFac, src1, src2, fac1, fac2)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: src1, src2
    real,    intent(in) :: fac1, fac2

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: src1Ptr, src2Ptr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr,  dst)
        call ml_memGetDataPtr(tileDesc, src1Ptr, src1)
        call ml_memGetDataPtr(tileDesc, src2Ptr, src2)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
            do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) =  dstFac*dstPtr(dstVars,i,j,k) &
                                            + fac1*src1Ptr(      :,i,j,k) &
                                            + fac2*src2Ptr(      :,i,j,k)
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, src2Ptr, src2)
        call ml_memReleaseDataPtr(tileDesc, src1Ptr, src1)
        call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars2

subroutine ml_memAddToVars3(dst, dstFac, src1, src2, src3, &
                             fac1, fac2, fac3)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: src1, src2, src3
    real,    intent(in) :: fac1, fac2, fac3

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: src1Ptr, src2Ptr, src3Ptr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr,  dst)
        call ml_memGetDataPtr(tileDesc, src1Ptr, src1)
        call ml_memGetDataPtr(tileDesc, src2Ptr, src2)
        call ml_memGetDataPtr(tileDesc, src3Ptr, src3)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
            do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) =   dstFac*dstPtr(dstVars,i,j,k) &
                                            + fac1*src1Ptr(      :,i,j,k) &
                                            + fac2*src2Ptr(      :,i,j,k) &
                                            + fac3*src3Ptr(      :,i,j,k)
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, src3Ptr, src3)
        call ml_memReleaseDataPtr(tileDesc, src2Ptr, src2)
        call ml_memReleaseDataPtr(tileDesc, src1Ptr, src1)
        call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars3

subroutine ml_memAddToVars4(dst, dstFac, src1, src2, src3, src4, &
                             fac1, fac2, fac3, fac4)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: src1, src2, src3, src4
    real,    intent(in) :: fac1, fac2, fac3, fac4

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: src1Ptr, src2Ptr, src3Ptr, src4Ptr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr,  dst)
        call ml_memGetDataPtr(tileDesc, src1Ptr, src1)
        call ml_memGetDataPtr(tileDesc, src2Ptr, src2)
        call ml_memGetDataPtr(tileDesc, src3Ptr, src3)
        call ml_memGetDataPtr(tileDesc, src4Ptr, src4)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
            do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) =  dstFac*dstPtr(dstVars,i,j,k) &
                                            + fac1*src1Ptr(      :,i,j,k) &
                                            + fac2*src2Ptr(      :,i,j,k) &
                                            + fac3*src3Ptr(      :,i,j,k) &
                                            + fac4*src4Ptr(      :,i,j,k)
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, src4Ptr, src4)
        call ml_memReleaseDataPtr(tileDesc, src3Ptr, src3)
        call ml_memReleaseDataPtr(tileDesc, src2Ptr, src2)
        call ml_memReleaseDataPtr(tileDesc, src1Ptr, src1)
        call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars4

subroutine ml_memAddToVars5(dst, dstFac, src1, src2, src3, src4, src5, &
                             fac1, fac2, fac3, fac4, fac5)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: src1, src2, src3, src4, src5
    real,    intent(in) :: fac1, fac2, fac3, fac4, fac5

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: src1Ptr, src2Ptr, src3Ptr, src4Ptr, src5Ptr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr,  dst)
        call ml_memGetDataPtr(tileDesc, src1Ptr, src1)
        call ml_memGetDataPtr(tileDesc, src2Ptr, src2)
        call ml_memGetDataPtr(tileDesc, src3Ptr, src3)
        call ml_memGetDataPtr(tileDesc, src4Ptr, src4)
        call ml_memGetDataPtr(tileDesc, src5Ptr, src5)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
            do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) =  dstFac*dstPtr(dstVars,i,j,k) &
                                            + fac1*src1Ptr(      :,i,j,k) &
                                            + fac2*src2Ptr(      :,i,j,k) &
                                            + fac3*src3Ptr(      :,i,j,k) &
                                            + fac4*src4Ptr(      :,i,j,k) &
                                            + fac5*src5Ptr(      :,i,j,k)
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, src5Ptr, src5)
        call ml_memReleaseDataPtr(tileDesc, src4Ptr, src4)
        call ml_memReleaseDataPtr(tileDesc, src3Ptr, src3)
        call ml_memReleaseDataPtr(tileDesc, src2Ptr, src2)
        call ml_memReleaseDataPtr(tileDesc, src1Ptr, src1)
        call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars5

subroutine ml_memAddToVars6(dst, dstFac, src1, src2, src3, src4, src5, src6, &
                             fac1, fac2, fac3, fac4, fac5, fac6)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: src1, src2, src3, src4, src5, src6
    real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                         src4Ptr, src5Ptr, src6Ptr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr,  dst)
        call ml_memGetDataPtr(tileDesc, src1Ptr, src1)
        call ml_memGetDataPtr(tileDesc, src2Ptr, src2)
        call ml_memGetDataPtr(tileDesc, src3Ptr, src3)
        call ml_memGetDataPtr(tileDesc, src4Ptr, src4)
        call ml_memGetDataPtr(tileDesc, src5Ptr, src5)
        call ml_memGetDataPtr(tileDesc, src6Ptr, src6)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
            do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) =  dstFac*dstPtr(dstVars,i,j,k) &
                                            + fac1*src1Ptr(      :,i,j,k) &
                                            + fac2*src2Ptr(      :,i,j,k) &
                                            + fac3*src3Ptr(      :,i,j,k) &
                                            + fac4*src4Ptr(      :,i,j,k) &
                                            + fac5*src5Ptr(      :,i,j,k) &
                                            + fac6*src6Ptr(      :,i,j,k)
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, src6Ptr, src6)
        call ml_memReleaseDataPtr(tileDesc, src5Ptr, src5)
        call ml_memReleaseDataPtr(tileDesc, src4Ptr, src4)
        call ml_memReleaseDataPtr(tileDesc, src3Ptr, src3)
        call ml_memReleaseDataPtr(tileDesc, src2Ptr, src2)
        call ml_memReleaseDataPtr(tileDesc, src1Ptr, src1)
        call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars6

subroutine ml_memAddToVars7(dst, dstFac, &
                             src1, src2, src3, src4, src5, src6, src7, &
                             fac1, fac2, fac3, fac4, fac5, fac6, fac7)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: src1, src2, src3, src4, src5, src6, src7
    real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                         src4Ptr, src5Ptr, src6Ptr, &
                                         src7Ptr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
    dstVars = MoL_unk_mask
    else
    dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr,  dst)
        call ml_memGetDataPtr(tileDesc, src1Ptr, src1)
        call ml_memGetDataPtr(tileDesc, src2Ptr, src2)
        call ml_memGetDataPtr(tileDesc, src3Ptr, src3)
        call ml_memGetDataPtr(tileDesc, src4Ptr, src4)
        call ml_memGetDataPtr(tileDesc, src5Ptr, src5)
        call ml_memGetDataPtr(tileDesc, src6Ptr, src6)
        call ml_memGetDataPtr(tileDesc, src7Ptr, src7)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
                do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) =  dstFac*dstPtr(dstVars,i,j,k) &
                                            + fac1*src1Ptr(      :,i,j,k) &
                                            + fac2*src2Ptr(      :,i,j,k) &
                                            + fac3*src3Ptr(      :,i,j,k) &
                                            + fac4*src4Ptr(      :,i,j,k) &
                                            + fac5*src5Ptr(      :,i,j,k) &
                                            + fac6*src6Ptr(      :,i,j,k) &
                                            + fac7*src7Ptr(      :,i,j,k)
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, src7Ptr, src7)
        call ml_memReleaseDataPtr(tileDesc, src6Ptr, src6)
        call ml_memReleaseDataPtr(tileDesc, src5Ptr, src5)
        call ml_memReleaseDataPtr(tileDesc, src4Ptr, src4)
        call ml_memReleaseDataPtr(tileDesc, src3Ptr, src3)
        call ml_memReleaseDataPtr(tileDesc, src2Ptr, src2)
        call ml_memReleaseDataPtr(tileDesc, src1Ptr, src1)
        call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars7

subroutine ml_memAddToVars8(dst, dstFac, &
                             src1, src2, src3, src4, src5, src6, src7, src8, &
                             fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: src1, src2, src3, src4, src5, src6, src7, src8
    real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                         src4Ptr, src5Ptr, src6Ptr, &
                                         src7Ptr, src8Ptr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr,  dst)
        call ml_memGetDataPtr(tileDesc, src1Ptr, src1)
        call ml_memGetDataPtr(tileDesc, src2Ptr, src2)
        call ml_memGetDataPtr(tileDesc, src3Ptr, src3)
        call ml_memGetDataPtr(tileDesc, src4Ptr, src4)
        call ml_memGetDataPtr(tileDesc, src5Ptr, src5)
        call ml_memGetDataPtr(tileDesc, src6Ptr, src6)
        call ml_memGetDataPtr(tileDesc, src7Ptr, src7)
        call ml_memGetDataPtr(tileDesc, src8Ptr, src8)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
            do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) =  dstFac*dstPtr(dstVars,i,j,k) &
                                    + fac1*src1Ptr(      :,i,j,k) &
                                    + fac2*src2Ptr(      :,i,j,k) &
                                    + fac3*src3Ptr(      :,i,j,k) &
                                    + fac4*src4Ptr(      :,i,j,k) &
                                    + fac5*src5Ptr(      :,i,j,k) &
                                    + fac6*src6Ptr(      :,i,j,k) &
                                    + fac7*src7Ptr(      :,i,j,k) &
                                    + fac8*src8Ptr(      :,i,j,k)
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, src8Ptr, src8)
        call ml_memReleaseDataPtr(tileDesc, src7Ptr, src7)
        call ml_memReleaseDataPtr(tileDesc, src6Ptr, src6)
        call ml_memReleaseDataPtr(tileDesc, src5Ptr, src5)
        call ml_memReleaseDataPtr(tileDesc, src4Ptr, src4)
        call ml_memReleaseDataPtr(tileDesc, src3Ptr, src3)
        call ml_memReleaseDataPtr(tileDesc, src2Ptr, src2)
        call ml_memReleaseDataPtr(tileDesc, src1Ptr, src1)
        call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars8

subroutine ml_memAddToVars9(dst, dstFac, &
                             src1, src2, src3, src4, src5, src6, src7, src8, src9, &
                             fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: src1, src2, src3, src4, src5, src6, src7, src8, src9
    real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                         src4Ptr, src5Ptr, src6Ptr, &
                                         src7Ptr, src8Ptr, src9Ptr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr,  dst)
        call ml_memGetDataPtr(tileDesc, src1Ptr, src1)
        call ml_memGetDataPtr(tileDesc, src2Ptr, src2)
        call ml_memGetDataPtr(tileDesc, src3Ptr, src3)
        call ml_memGetDataPtr(tileDesc, src4Ptr, src4)
        call ml_memGetDataPtr(tileDesc, src5Ptr, src5)
        call ml_memGetDataPtr(tileDesc, src6Ptr, src6)
        call ml_memGetDataPtr(tileDesc, src7Ptr, src7)
        call ml_memGetDataPtr(tileDesc, src8Ptr, src8)
        call ml_memGetDataPtr(tileDesc, src9Ptr, src9)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
            do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) =  dstFac*dstPtr(dstVars,i,j,k) &
                                    + fac1*src1Ptr(      :,i,j,k) &
                                    + fac2*src2Ptr(      :,i,j,k) &
                                    + fac3*src3Ptr(      :,i,j,k) &
                                    + fac4*src4Ptr(      :,i,j,k) &
                                    + fac5*src5Ptr(      :,i,j,k) &
                                    + fac6*src6Ptr(      :,i,j,k) &
                                    + fac7*src7Ptr(      :,i,j,k) &
                                    + fac8*src8Ptr(      :,i,j,k) &
                                    + fac9*src9Ptr(      :,i,j,k)
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, src9Ptr, src9)
        call ml_memReleaseDataPtr(tileDesc, src8Ptr, src8)
        call ml_memReleaseDataPtr(tileDesc, src7Ptr, src7)
        call ml_memReleaseDataPtr(tileDesc, src6Ptr, src6)
        call ml_memReleaseDataPtr(tileDesc, src5Ptr, src5)
        call ml_memReleaseDataPtr(tileDesc, src4Ptr, src4)
        call ml_memReleaseDataPtr(tileDesc, src3Ptr, src3)
        call ml_memReleaseDataPtr(tileDesc, src2Ptr, src2)
        call ml_memReleaseDataPtr(tileDesc, src1Ptr, src1)
        call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars9

subroutine ml_memAddToVars10(dst, dstFac, &
                              src1, src2, src3, src4, src5, src6, src7, src8, &
                              src9, src10, &
                              fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, &
                              fac9, fac10)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: src1, src2, src3, src4, src5, src6, src7, src8, src9, src10
    real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9, fac10

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                         src4Ptr, src5Ptr, src6Ptr, &
                                         src7Ptr, src8Ptr, src9Ptr, &
                                         src10Ptr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr,  dst)
        call ml_memGetDataPtr(tileDesc, src1Ptr, src1)
        call ml_memGetDataPtr(tileDesc, src2Ptr, src2)
        call ml_memGetDataPtr(tileDesc, src3Ptr, src3)
        call ml_memGetDataPtr(tileDesc, src4Ptr, src4)
        call ml_memGetDataPtr(tileDesc, src5Ptr, src5)
        call ml_memGetDataPtr(tileDesc, src6Ptr, src6)
        call ml_memGetDataPtr(tileDesc, src7Ptr, src7)
        call ml_memGetDataPtr(tileDesc, src8Ptr, src8)
        call ml_memGetDataPtr(tileDesc, src9Ptr, src9)
        call ml_memGetDataPtr(tileDesc, src10Ptr, src10)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
            do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) =  dstFac*dstPtr(dstVars,i,j,k) &
                                    + fac1*src1Ptr(      :,i,j,k) &
                                    + fac2*src2Ptr(      :,i,j,k) &
                                    + fac3*src3Ptr(      :,i,j,k) &
                                    + fac4*src4Ptr(      :,i,j,k) &
                                    + fac5*src5Ptr(      :,i,j,k) &
                                    + fac6*src6Ptr(      :,i,j,k) &
                                    + fac7*src7Ptr(      :,i,j,k) &
                                    + fac8*src8Ptr(      :,i,j,k) &
                                    + fac9*src9Ptr(      :,i,j,k) &
                                    + fac10*src10Ptr(    :,i,j,k)
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, src10Ptr, src10)
        call ml_memReleaseDataPtr(tileDesc, src9Ptr, src9)
        call ml_memReleaseDataPtr(tileDesc, src8Ptr, src8)
        call ml_memReleaseDataPtr(tileDesc, src7Ptr, src7)
        call ml_memReleaseDataPtr(tileDesc, src6Ptr, src6)
        call ml_memReleaseDataPtr(tileDesc, src5Ptr, src5)
        call ml_memReleaseDataPtr(tileDesc, src4Ptr, src4)
        call ml_memReleaseDataPtr(tileDesc, src3Ptr, src3)
        call ml_memReleaseDataPtr(tileDesc, src2Ptr, src2)
        call ml_memReleaseDataPtr(tileDesc, src1Ptr, src1)
        call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars10

subroutine ml_memAddToVars11(dst, dstFac, &
                              src1, src2, src3, src4, src5, src6, src7, src8, &
                              src9, src10, src11, &
                              fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, &
                              fac9, fac10, fac11)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                           src7, src8, src9, src10, src11
    real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                           fac7, fac8, fac9, fac10, fac11

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                         src4Ptr, src5Ptr, src6Ptr, &
                                         src7Ptr, src8Ptr, src9Ptr, &
                                         src10Ptr, src11Ptr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr,  dst)
        call ml_memGetDataPtr(tileDesc, src1Ptr, src1)
        call ml_memGetDataPtr(tileDesc, src2Ptr, src2)
        call ml_memGetDataPtr(tileDesc, src3Ptr, src3)
        call ml_memGetDataPtr(tileDesc, src4Ptr, src4)
        call ml_memGetDataPtr(tileDesc, src5Ptr, src5)
        call ml_memGetDataPtr(tileDesc, src6Ptr, src6)
        call ml_memGetDataPtr(tileDesc, src7Ptr, src7)
        call ml_memGetDataPtr(tileDesc, src8Ptr, src8)
        call ml_memGetDataPtr(tileDesc, src9Ptr, src9)
        call ml_memGetDataPtr(tileDesc, src10Ptr, src10)
        call ml_memGetDataPtr(tileDesc, src11Ptr, src11)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
            do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) =  dstFac*dstPtr(dstVars,i,j,k) &
                                    + fac1*src1Ptr(      :,i,j,k) &
                                    + fac2*src2Ptr(      :,i,j,k) &
                                    + fac3*src3Ptr(      :,i,j,k) &
                                    + fac4*src4Ptr(      :,i,j,k) &
                                    + fac5*src5Ptr(      :,i,j,k) &
                                    + fac6*src6Ptr(      :,i,j,k) &
                                    + fac7*src7Ptr(      :,i,j,k) &
                                    + fac8*src8Ptr(      :,i,j,k) &
                                    + fac9*src9Ptr(      :,i,j,k) &
                                    + fac10*src10Ptr(    :,i,j,k) &
                                    + fac11*src11Ptr(    :,i,j,k)
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, src11Ptr, src11)
        call ml_memReleaseDataPtr(tileDesc, src10Ptr, src10)
        call ml_memReleaseDataPtr(tileDesc, src9Ptr, src9)
        call ml_memReleaseDataPtr(tileDesc, src8Ptr, src8)
        call ml_memReleaseDataPtr(tileDesc, src7Ptr, src7)
        call ml_memReleaseDataPtr(tileDesc, src6Ptr, src6)
        call ml_memReleaseDataPtr(tileDesc, src5Ptr, src5)
        call ml_memReleaseDataPtr(tileDesc, src4Ptr, src4)
        call ml_memReleaseDataPtr(tileDesc, src3Ptr, src3)
        call ml_memReleaseDataPtr(tileDesc, src2Ptr, src2)
        call ml_memReleaseDataPtr(tileDesc, src1Ptr, src1)
        call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars11

subroutine ml_memAddToVars12(dst, dstFac, &
                              src1, src2, src3, src4, src5, src6, src7, src8, &
                              src9, src10, src11, src12, &
                              fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, &
                              fac9, fac10, fac11, fac12)
    use MoL_variables,   only: MoL_nvars, MoL_unk_mask, MoL_scratch_mask
    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr

    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_tile,      only: Grid_tile_t
    use Grid_iterator,  only: Grid_iterator_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: dst
    real,    intent(in) :: dstFac
    integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                           src7, src8, src9, src10, src11, src12
    real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                           fac7, fac8, fac9, fac10, fac11, fac12

    integer, dimension(MoL_nvars) :: dstVars
    real, dimension(:,:,:,:), pointer :: dstPtr
    real, dimension(:,:,:,:), pointer :: src1Ptr, src2Ptr, src3Ptr, &
                                         src4Ptr, src5Ptr, src6Ptr, &
                                         src7Ptr, src8Ptr, src9Ptr, &
                                         src10Ptr, src11Ptr, src12Ptr

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    integer :: i, j, k

    if (dst .eq. MOL_EVOLVED) then
        dstVars = MoL_unk_mask
    else
        dstVars = MoL_scratch_mask
    end if

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, dstPtr,  dst)
        call ml_memGetDataPtr(tileDesc, src1Ptr, src1)
        call ml_memGetDataPtr(tileDesc, src2Ptr, src2)
        call ml_memGetDataPtr(tileDesc, src3Ptr, src3)
        call ml_memGetDataPtr(tileDesc, src4Ptr, src4)
        call ml_memGetDataPtr(tileDesc, src5Ptr, src5)
        call ml_memGetDataPtr(tileDesc, src6Ptr, src6)
        call ml_memGetDataPtr(tileDesc, src7Ptr, src7)
        call ml_memGetDataPtr(tileDesc, src8Ptr, src8)
        call ml_memGetDataPtr(tileDesc, src9Ptr, src9)
        call ml_memGetDataPtr(tileDesc, src10Ptr, src10)
        call ml_memGetDataPtr(tileDesc, src11Ptr, src11)
        call ml_memGetDataPtr(tileDesc, src12Ptr, src12)

        do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
            do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
                do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                    dstPtr(dstVars,i,j,k) =  dstFac*dstPtr(dstVars,i,j,k) &
                                    + fac1*src1Ptr(      :,i,j,k) &
                                    + fac2*src2Ptr(      :,i,j,k) &
                                    + fac3*src3Ptr(      :,i,j,k) &
                                    + fac4*src4Ptr(      :,i,j,k) &
                                    + fac5*src5Ptr(      :,i,j,k) &
                                    + fac6*src6Ptr(      :,i,j,k) &
                                    + fac7*src7Ptr(      :,i,j,k) &
                                    + fac8*src8Ptr(      :,i,j,k) &
                                    + fac9*src9Ptr(      :,i,j,k) &
                                    + fac10*src10Ptr(    :,i,j,k) &
                                    + fac11*src11Ptr(    :,i,j,k) &
                                    + fac12*src12Ptr(    :,i,j,k)
                end do ! i
            end do ! j
        end do ! k

        call ml_memReleaseDataPtr(tileDesc, src12Ptr, src12)
        call ml_memReleaseDataPtr(tileDesc, src11Ptr, src11)
        call ml_memReleaseDataPtr(tileDesc, src10Ptr, src10)
        call ml_memReleaseDataPtr(tileDesc, src9Ptr, src9)
        call ml_memReleaseDataPtr(tileDesc, src8Ptr, src8)
        call ml_memReleaseDataPtr(tileDesc, src7Ptr, src7)
        call ml_memReleaseDataPtr(tileDesc, src6Ptr, src6)
        call ml_memReleaseDataPtr(tileDesc, src5Ptr, src5)
        call ml_memReleaseDataPtr(tileDesc, src4Ptr, src4)
        call ml_memReleaseDataPtr(tileDesc, src3Ptr, src3)
        call ml_memReleaseDataPtr(tileDesc, src2Ptr, src2)
        call ml_memReleaseDataPtr(tileDesc, src1Ptr, src1)
        call ml_memReleaseDataPtr(tileDesc, dstPtr,  dst)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_memAddToVars12
