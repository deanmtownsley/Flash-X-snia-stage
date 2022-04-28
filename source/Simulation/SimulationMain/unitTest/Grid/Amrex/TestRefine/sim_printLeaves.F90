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

subroutine sim_printLeaves(title)
    use gr_amrexInterface, ONLY : gr_getFinestLevel
    use Simulation_data,   ONLY : leaves, &
                                  MIN_REFINE_LEVEL

    implicit none

    character(*), intent(IN) :: title

    integer :: level, j
    integer :: finest_level

    write(*,'(A)') title
    write(*,'(A)') "-----------------------------------"
    call gr_getFinestLevel(finest_level)
    do level = MIN_REFINE_LEVEL, finest_level 
        if (.NOT. allocated(leaves(level)%blocks))   CYCLE

        associate(blks => leaves(level)%blocks)
            write(*,'(A,I2)') "Leaf blocks at level ", level
            do j = 1, SIZE(blks, 1)
                write(*,'(A,I3,A,I3,A,I3,A,I3,A)') "     From (", &
                         blks(j, 1), ", ", &
                         blks(j, 2), ") to (", &
                         blks(j, 3), ", ", &
                         blks(j, 4), ")"
            end do
        end associate

    end do
end subroutine sim_printLeaves

