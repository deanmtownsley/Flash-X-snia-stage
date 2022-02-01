!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestInit/Simulation_initBlock
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
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer (IN) ::blockId, 
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes the Grid with a composit number which is a combination
!!  of the block number and the indices of the cell
!! 
!! ARGUMENTS
!!
!!  blockId -          the blockId to update
!!  
!!
!!
!!***

#include "Simulation.h"
#include "constants.h"
#include "sim_constants.h"

subroutine Simulation_initBlock(initData, tileDesc)
    use Grid_tile,     ONLY : Grid_tile_t 
    use sim_interface, ONLY : sim_writeDataPoints

    implicit none

    real,                         pointer :: initData(:, :, :, :)
    type(Grid_tile_t), intent(IN)         :: tileDesc

    real    :: points(2, 2)
    real    :: values(2)

    integer :: i, j, k, var

    associate(lo => tileDesc%limits(LOW,  :), &
              hi => tileDesc%limits(HIGH, :))
        do           var = UNK_VARS_BEGIN, UNK_VARS_END
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        initData(i, j, k, var) = 0.0
                    end do
                end do
            end do
        end do
    end associate

    points(:, :) = 0.0
    points(1, :) = [0.16, 0.67]
    points(2, :) = [0.11, 0.38]
    values(:) = 0.0
    values(1) = REFINE_TO_L3
    values(2) = REFINE_TO_L2

    call sim_writeDataPoints(initData, tileDesc, points, values)
end subroutine Simulation_initBlock

