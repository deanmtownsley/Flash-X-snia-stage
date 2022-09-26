!> @file
!! @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock

#include "constants.h"

!> @ingroup GridMilhoja
!!
!! @brief Ask Milhoja to prepare data needed by IO unit for writing to file
!!
!! @details
!! A routine that should be called in preparation for checkpointing so that
!! Milhoja can prepare data that the IO unit needs for writing data to file.
!!
!! In particular, this updates some arrays that provide a rudimentary simulacrum
!! of same data private to the PARAMESH implementation.  It should be called
!! immediately after Grid finishes refining or updating refinement.
!!
!! @attention
!! This is a quick & dirty implementation built up from the AMReX version.
!! In particular, it will only work for pseudo-UG runs.  This is acceptable
!! as Milhoja only works in pseudo-UG presently.
!!
!! @todo Finish implementation and update for AMR when this is available
!! through Milhoja.
subroutine gr_updateDataForIo()
    use Driver_interface, ONLY : Driver_abort
    use Grid_interface,   ONLY : Grid_getTileIterator, &
                                 Grid_releaseTileIterator, &
                                 Grid_getLocalNumBlks
    use Grid_iterator,    ONLY : Grid_iterator_t
    use Grid_tile,        ONLY : Grid_tile_t
    use gr_specificData,  ONLY : gr_ioLocalNumBlocks, &
                                 gr_ioGlobalNumBlocks, &
                                 gr_ioBlkNodeType, &
                                 gr_ioBlkCoords, &
                                 gr_ioBlkBsize , &
                                 gr_ioBlkBoundBox, &
                                 gr_ioBlkLrefine

    implicit none

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t)     :: tileDesc
    integer               :: nBlocks
    integer               :: j

    CALL Grid_getLocalNumBlks(nBlocks)
    gr_ioLocalNumBlocks = nBlocks

    if (ALLOCATED(gr_ioBlkLrefine)) then
        deallocate(gr_ioBlkLrefine)
    end if
    if (ALLOCATED(gr_ioBlkNodeType)) then
        deallocate(gr_ioBlkNodeType)
    end if
    if (ALLOCATED(gr_ioBlkCoords)) then
        deallocate(gr_ioBlkCoords)
    end if
    if (ALLOCATED(gr_ioBlkBsize)) then
        deallocate(gr_ioBlkBsize)
    end if
    if (ALLOCATED(gr_ioBlkBoundBox)) then
        deallocate(gr_ioBlkBoundBox)
    end if

    allocate(gr_ioBlkLrefine(nBlocks))
    allocate(gr_ioBlkNodeType(nBlocks))
    allocate(gr_ioBlkCoords(            MDIM, nBlocks))
    allocate(gr_ioBlkBsize(             MDIM, nBlocks))
    allocate(gr_ioBlkBoundBox(LOW:HIGH, MDIM, nBlocks))

    j = 1
    CALL Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.)
    do while (itor%isValid())
        CALL itor%currentTile(tileDesc)

        CALL tileDesc%boundBox(gr_ioBlkBoundBox(:, :, j))

        gr_ioBlkNodeType(j)  = LEAF
        gr_ioBlkLrefine(j)   = tileDesc%level
        gr_ioBlkBsize(:, j)  =        gr_ioBlkBoundBox(HIGH, :, j) - gr_ioBlkBoundBox(LOW, :, j)
        gr_ioBlkCoords(:, j) = 0.5 * (gr_ioBlkBoundBox(HIGH, :, j) + gr_ioBlkBoundBox(LOW, :, j))

        j = j + 1
        CALL itor%next()
    enddo
    CALL Grid_releaseTileIterator(itor)

    if (gr_ioLocalNumBlocks .NE. (j - 1)) then
        call Driver_abort("[gr_updateDataForIo] Inconsistent N blocks")
    end if

    ! This is not set in this routine, so give it a garbage default value.
    gr_ioGlobalNumBlocks = -HUGE(gr_ioGlobalNumBlocks)
end subroutine gr_updateDataForIo

