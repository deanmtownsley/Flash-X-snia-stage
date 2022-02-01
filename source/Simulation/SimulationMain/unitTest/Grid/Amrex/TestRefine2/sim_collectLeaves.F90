#include "constants.h"
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

subroutine sim_collectLeaves
    use gr_amrexInterface,    ONLY : gr_getFinestLevel
    use block_iterator,       ONLY : block_iterator_t 
    use Grid_tile,       ONLY : Grid_tile_t 
    use Simulation_data,      ONLY : leaves, &
                                     MIN_REFINE_LEVEL, MAX_REFINE_LEVEL

    type(block_iterator_t) :: itor
    type(Grid_tile_t) :: blockDesc

    logical :: gridChanged
    integer :: finest_level
    integer :: block_count
    integer :: lev, j

    ! Regenerate leaf block data structure
    call gr_getFinestLevel(finest_level)
    do lev = MIN_REFINE_LEVEL, MAX_REFINE_LEVEL
        block_count = 0
        itor = block_iterator_t(LEAF, level=lev)
        do while (itor%is_valid())
            call itor%blkMetaData(blockDesc)
 
            block_count = block_count + 1

            call itor%next()
        end do

        if (allocated(leaves(lev)%blocks)) then
            deallocate(leaves(lev)%blocks)
        end if

        if (block_count > 0) then
            allocate(leaves(lev)%blocks(block_count, 4))
        end if
    end do

    ! Populate leaf block data structure
    do lev = MIN_REFINE_LEVEL, finest_level 
        itor = block_iterator_t(LEAF, level=lev)

        j = 1
        do while (itor%is_valid())
            call itor%blkMetaData(blockDesc)

            associate(lo => blockDesc%limits(LOW,  :), &
                      hi => blockDesc%limits(HIGH, :))
                leaves(lev)%blocks(j, :) = [lo(IAXIS), lo(JAXIS), &
                                            hi(IAXIS), hi(JAXIS)]
            end associate

            j = j + 1
            call itor%next()
        end do
    end do
end subroutine sim_collectLeaves

