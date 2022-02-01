subroutine sim_advance(step, points, values, set_msg, leaf_msg)
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
    use Grid_interface,       ONLY : Grid_updateRefinement, &
                                     Grid_getBlkPtr, Grid_releaseBlkPtr
    use gr_amrexInterface,    ONLY : gr_getFinestLevel, &
                                     gr_restrictAllLevels
    use block_iterator,       ONLY : block_iterator_t 
    use Grid_tile,       ONLY : Grid_tile_t 
    use Driver_interface,     ONLY : Driver_abort
    use sim_interface,        ONLY : sim_writeDataPoints, &
                                     sim_collectLeaves, &
                                     sim_printLeaves

    implicit none

#include "constants.h"
 
    integer,      intent(IN)    :: step
    real,         intent(IN)    :: points(:, :)
    real,         intent(IN)    :: values(:)
    character(*), intent(IN)    :: set_msg
    character(*), intent(IN)    :: leaf_msg

    real, contiguous, pointer :: solnData(:,:,:,:)
 
    type(block_iterator_t) :: itor
    type(Grid_tile_t) :: blockDesc

    logical :: gridChanged
    integer :: finest_level
    integer :: lev

    !!!!! ZERO DATA AND WRITE GIVEN POINTS
    write(*,*)
    write(*,'(A,I2,A,I2)') "ADVANCE STEPS ", step, "/", step+1
    write(*,'(A)') "--------------------------------------------------------------"
    write(*,*) set_msg
    write(*,*)

    ! Write to leaf blocks first.  AMReX level indexing is 0-based
    call gr_getFinestLevel(finest_level)
    do lev = 1, finest_level
        itor = block_iterator_t(LEAF, level=lev)
        do while (itor%is_valid())
            call itor%blkMetaData(blockDesc)
            call Grid_getBlkPtr(blockDesc, solnData)

            solnData = 0.0d0
            call sim_writeDataPoints(solnData, blockDesc,  points, values)

            call Grid_releaseBlkPtr(blockDesc, solnData)
            call itor%next()
        end do
    end do

    ! Propogate leaf data to coarse, non-leaf blocks
    call gr_restrictAllLevels(CENTER, convertPtoC=.FALSE., convertCtoP=.FALSE.)

    !!!!! REFINE/DEREFINE BASED ON NEW DATA
    write(*,*)
    write(*,*) "EXECUTING REFINE/DEREFINE"
    write(*,*)

    ! Should only refine on every other step (nrefs = 2)
    gridChanged = .FALSE.
    call Grid_updateRefinement(step,   DBLE(step  ), gridChanged) 
    if (gridChanged) then
        call Driver_abort("[sim_advance] Should not refine on odd steps")
    end if

    call Grid_updateRefinement(step+1, DBLE(step+1), gridChanged) 
    if (.NOT. gridChanged) then
        call Driver_abort("[sim_advance] Should refine on even steps")
    end if

    call sim_collectLeaves
end subroutine sim_advance

