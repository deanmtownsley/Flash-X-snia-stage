!!****ih* source/Grid/GridMain/AMR/Paramesh4/Grid_iterator
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
!!
!!
!!****

!! defines IMPURE_ELEMENTAL:
#include "FortranLangFeatures.fh"
#include "Simulation.h"
#include "constants.h"

module Grid_iterator
    use tree, ONLY : lnblocks

    implicit none

    private
    public :: build_iterator, destroy_iterator

    !!****ic* Grid_iterator/Grid_iterator_t
    !!
    !! NAME
    !!  Grid_iterator
    !!
    !!****
    type, public :: Grid_iterator_t
        integer :: curBlk   = 0
        integer :: nodetype = LEAF 
        integer :: lev      = INVALID_LEVEL
        integer, private :: tileSize(MDIM) = (/1,1,1/)
        integer, private :: nxt=1,nyt=1,nzt=1
        logical, private :: tiling         = .FALSE.
        integer, private :: curTile        = 0 ! local tileNo in the block, 0-based
        integer, private :: ntiles         = 0 ! number of local tiles per block
        integer, private :: nthreads       = 1 ! number of threads to participate in iterations
#ifdef _OPENMP
        integer, private :: threadOffs     = 0
#endif
    contains
        procedure, public :: isValid
        procedure         :: next_ ! simple next
        procedure         :: nextN ! next by n
        procedure, public :: next
        procedure, public :: currentTile
    end type Grid_iterator_t

contains

    !!****im* Grid_iterator/build_iterator
    !!
    !! NAME
    !!  build_iterator
    !!
    !! SYNOPOSIS
    !!  build_iterator(Grid_iterator_t(OUT)    :: itor,
    !!                 integer(IN)           :: nodetype,
    !!                 integer(IN), optional :: level,
    !!                 logical(IN), optional :: tiling,
    !!                 integer(IN), optional :: nthreads)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of blocks or
    !!  tiles within the current paramesh octree structure.  The iterator is already
    !!  set to the first matching block/tile.
    !!  
    !! ARGUMENTS
    !!  itor     - the iterator object
    !!  nodetype - the class of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS).
    !!             Refer to the documentation for the Paramesh version of
    !!             Grid_getTileIterator for more details.
    !!  level    - iterate only over all blocks/tiles of the correct nodetype
    !!             that are located at this level of refinement.  Refer to the
    !!             documentation for the Paramesh version of Grid_getTileIterator
    !!             for more details.  A level value of UNSPEC_LEVEL is equivalent
    !!             to omitting this optional argument.
    !!  tiling   - an optional optimization hint.  If TRUE, proper tiling is
    !!             requested.
    !!  nthreads - an optional argument that may affect division of the iteration
    !!             space in an active OpenMP parallel region.
    !!             Can usually be omitted since the default behavior is reasonable.
    !!
    !! SEE ALSO
    !!  constants.h
    !!****
    subroutine build_iterator(itor, nodetype, level, tiling, tileSize, nthreads)
        use Driver_interface, ONLY : Driver_abort
        !$ use omp_lib,       ONLY : omp_get_num_threads, omp_get_thread_num
        type(Grid_iterator_t), intent(OUT) :: itor
        integer,               intent(IN)  :: nodetype
        integer,               intent(IN)  :: level
        logical,               intent(IN)  :: tiling
        integer,               intent(IN)  :: tileSize(1:MDIM)
        integer,               intent(IN), optional :: nthreads

        if (present(nthreads)) then
           itor%nthreads = nthreads
        else
#ifdef _OPENMP
           itor%nthreads = omp_get_num_threads()
#else
           itor%nthreads = 1
#endif
        end if
#ifdef _OPENMP
        itor%threadOffs = omp_get_thread_num()
        if (itor%threadOffs .GE. itor%nthreads) itor%threadOffs = 0
#endif

        itor%nodetype = nodetype
        itor%lev = level
        itor%tiling = tiling
        itor%tileSize(1:NDIM) = tileSize(1:NDIM)
        if (tiling) then
           itor%nxt = max(1,(NXB+tileSize(IAXIS)-1)/tileSize(IAXIS))
           itor%ntiles = itor%nxt
#if NDIM > 1
           itor%nyt = max(1,(NYB+tileSize(JAXIS)-1)/tileSize(JAXIS))
           itor%ntiles = itor%ntiles * itor%nyt
#endif
#if NDIM >= 3
           itor%nzt = max(1,(NZB+tileSize(KAXIS)-1)/tileSize(KAXIS))
           itor%ntiles = itor%ntiles * itor%nzt
#endif
        else
          itor% ntiles = 1
        end if
        itor%curBlk = 0
        itor%curTile = itor%ntiles - 1
#ifdef _OPENMP
        call itor%nextN(1 + itor%threadOffs)
#else
        call itor%next_()
#endif
    end subroutine build_iterator

    !!****im* Grid_iterator_t/destroy_iterator
    !!
    !! NAME
    !!  destroy_iterator
    !!
    !! SYNPOSIS
    !!  Destroy given iterator.
    !!
    !! DESCRIPTION
    !!  Clean-up block interator object at destruction
    !!
    !!****
    IMPURE_ELEMENTAL subroutine destroy_iterator(itor)
        type(Grid_iterator_t), intent(INOUT) :: itor

        itor%curBlk = lnblocks + 1
    end subroutine destroy_iterator

    !!****m* Grid_iterator_t/isValid
    !!
    !! NAME
    !!  isValid
    !!
    !! SYNPOSIS
    !!  logical valid = itor%isValid()
    !!
    !! DESCRIPTION
    !!  Determine if the iterator is currently set to a valid block.
    !!
    !! RETURN VALUE 
    !!  True if iterator is currently set to a valid block
    !!
    !!****
    logical function isValid(this)
        class(Grid_iterator_t), intent(IN) :: this

        isValid = (this%curBlk <= lnblocks)
    end function isValid

    !!****m* Grid_iterator_t/next_
    !!
    !! NAME
    !!  next_
    !!
    !! SYNPOSIS
    !!  call itor%next_()
    !!
    !! DESCRIPTION
    !!  Advance the iterator to the next tile managed by the process that meets
    !!  the iterator constraints given at instantiation.
    !!
    !!****
    subroutine next_(this)
        use gr_parameshInterface, ONLY : gr_blockMatch

        class(Grid_iterator_t), intent(INOUT) :: this

        integer :: j

        if (this%curTile < this%ntiles-1) then
           this%curTile = this%curTile + 1
        else
           if (this%lev == UNSPEC_LEVEL) then
            ! No level given at creation
              do j = this%curBlk + 1, lnblocks
                 if (gr_blockMatch(j, this%nodetype)) EXIT
              end do
           else
              do j = this%curBlk + 1, lnblocks
                 if (gr_blockMatch(j, this%nodetype, this%lev)) EXIT
              end do
           end if
           this%curBlk = j
           this%curTile = 0  
        end if

    end subroutine next_

    !!****m* Grid_iterator_t/nextN
    !!
    !! NAME
    !!  nextN
    !!
    !! SYNPOSIS
    !!  call itor%nextN(n)
    !!
    !! DESCRIPTION
    !!  Advance the iterator up to n times, as if next() was invoked
    !!  until either n invocations were done or isValid becomes false,
    !!  whichever occurs first.
    !!
    !!****
    subroutine nextN(this, n)
        class(Grid_iterator_t), intent(INOUT) :: this
        integer,                intent(IN)    :: n

        integer :: i

        do i = 1,n
           if (this%curBlk > lnblocks) EXIT
           call this%next_()
        end do

    end subroutine nextN

    !!****m* Grid_iterator_t/next
    !!
    !! NAME
    !!  next
    !!
    !! SYNPOSIS
    !!  call itor%next()
    !!
    !! DESCRIPTION
    !!  Advance the iterator to the next tile managed by the current process and thread
    !! that meets the iterator constraints given at instantiation.
    !!
    !!****
    subroutine next(this)
        class(Grid_iterator_t), intent(INOUT) :: this

#ifdef _OPENMP
        call this%nextN(this%nthreads)
#else
        call this%next_()
#endif

    end subroutine next

    !!****m* Grid_iterator_t/blkMetaData
    !!
    !! NAME
    !!  blkMetaData 
    !!
    !! SYNPOSIS
    !!  call itor%currentTile(Grid_tile_t(OUT) : tileDesc)
    !!
    !! DESCRIPTION
    !!  Obtain meta data that characterizes the block currently set in the
    !!  iterator.
    !!
    !!****
    subroutine currentTile(this, tileDesc)
        use gr_specificData, ONLY : gr_oneBlock
        use Grid_tile,       ONLY : Grid_tile_t 
        use tree,            ONLY : lrefine, &
                                    lrefine_max

        class(Grid_iterator_t), intent(IN)  :: this
        type(Grid_tile_t),      intent(OUT) :: tileDesc

        integer :: cornerID(1:MDIM)
        integer :: loBnd(1:MDIM), hiBnd(1:MDIM)
        integer :: blkLim(LOW:HIGH, 1:MDIM)
        integer :: tx,ty,tz

        tileDesc%id = this%curBlk
        tileDesc%level = lrefine(tileDesc%id)

        tileDesc%cid = gr_oneBlock(tileDesc%id)%cornerID
        tileDesc%stride = 2**(lrefine_max - tileDesc%level)

        associate(lo      => tileDesc%limits(LOW, :), &
                  hi      => tileDesc%limits(HIGH, :), &
                  loGrown => tileDesc%grownLimits(LOW, :), &
                  hiGrown => tileDesc%grownLimits(HIGH, :), &
                  loBlkGC => tileDesc%blkLimitsGC(LOW, :), &
                  hiBlkGC => tileDesc%blkLimitsGC(HIGH, :), &
                  blkId   => tileDesc%id, &
                  cid     => tileDesc%cid)
            cornerID = (cid - 1) / 2**(lrefine_max-lrefine(blkID)) + 1

            ! Get limits with block-local indexing of cells
            blkLim(LOW,  IAXIS) = GRID_ILO
            blkLim(HIGH, IAXIS) = GRID_IHI
            blkLim(LOW,  JAXIS) = GRID_JLO
            blkLim(HIGH, JAXIS) = GRID_JHI
            blkLim(LOW,  KAXIS) = GRID_KLO
            blkLim(HIGH, KAXIS) = GRID_KHI

            ! Convert to global indexing of cells
            loBlkGC(:) = blkLim(LOW,  :) - 1 + cornerID(:)
            hiBlkGC(:) = blkLim(HIGH, :) - 1 + cornerID(:)
            loBlkGC(1:NDIM) = loBlkGC(1:NDIM) - (NGUARD+NGUARD)

            call tile2txtytz(this%curTile,tx,ty,tz)
            ! Convert to global indexing of cells
            lo(:)   = blkLim(LOW,  :) - 1 + cornerID(:)
            hi(:)   = blkLim(HIGH, :) - 1 + cornerID(:)
            lo(1:NDIM)   = lo(1:NDIM) - NGUARD
            hi(1:NDIM)   = hi(1:NDIM) - NGUARD

            hi(IAXIS) = min(hi(IAXIS), lo(IAXIS)+NXB*(tx+1)/this%nxt-1)
            lo(IAXIS) =     lo(IAXIS) + NXB*tx/this%nxt
#if (NDIM > 1)
            hi(JAXIS) = min(hi(JAXIS), lo(JAXIS)+NYB*(ty+1)/this%nyt-1)
            lo(JAXIS) =     lo(JAXIS) + NYB*ty/this%nyt
#endif
#if (NDIM > 2)
            hi(KAXIS) = min(hi(KAXIS), lo(KAXIS)+NZB*(tz+1)/this%nzt-1)
            lo(KAXIS) =     lo(KAXIS) + NZB*tz/this%nzt
#endif
            
            ! As long as there is no tiling with Paramesh, the grown tile
            ! is just the block + GC halo
            loGrown(:) = loBlkGC(:)
            hiGrown(:) = hiBlkGC(:)
            if (this%tiling) then
               loBnd(:) = max(0, 1 - (/tx,ty,tz/) ) ! indicates whether range includes first cell
!!$            hiBnd(:) = 1 - min(1, (/this%nxt-1,this%nyt-1,this%nzt-1/) - (/tx,ty,tz/) )
               hiBnd(:) = max(0, (/tx,ty,tz/) - (/this%nxt-2,this%nyt-2,this%nzt-2/) )
               loGrown(1:NDIM) = max(loBlkGC(1:NDIM), lo(1:NDIM) - loBnd(1:NDIM)*NGUARD)
               hiGrown(1:NDIM) = min(hiBlkGC(1:NDIM), hi(1:NDIM) + hiBnd(1:NDIM)*NGUARD)
            end if
        end associate

      contains
        subroutine tile2txtytz(tileID,tx,ty,tz)
          integer,intent(IN) :: tileID
          integer,intent(OUT) :: tx,ty,tz

          integer :: tn
          tn = tileID

          associate(nxt      => this%nxt, &
                    nyt      => this%nyt, &
                    nzt      => this%nzt)
            tx = mod(tn,nxt)
            ty = mod(tn/nxt,nyt)
            tz = tn/(nxt*nyt)
          end associate
        end subroutine tile2txtytz

    end subroutine currentTile

end module Grid_iterator

