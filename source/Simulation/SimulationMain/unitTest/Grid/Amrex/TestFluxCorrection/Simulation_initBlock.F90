!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestFluxCorrection/Simulation_initBlock
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
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!  Simulation_initBlock(       real :: initData,
!!                       integer(IN) :: tileDesc) 
!!
!! DESCRIPTION
!!  Initializes the data to the 1-based level at which the block exists
!! 
!! ARGUMENTS
!!  initData - the cell-centered data structure to which data is written
!!  tileDesc - index of block whose cell-centered data is to be initialized 
!!
!!***

#include "Simulation.h"
#include "constants.h"

subroutine Simulation_initBlock(initData, tileDesc)
    use Grid_tile, ONLY : Grid_tile_t 

    implicit none

    real,                         pointer :: initData(:, :, :, :)
    type(Grid_tile_t), intent(IN)         :: tileDesc

    integer :: i, j, k, var

    ! Initialize data.  The values are not important for this test
    associate(lo => tileDesc%limits(LOW,  :), &
              hi => tileDesc%limits(HIGH, :))
        do           var = UNK_VARS_BEGIN, UNK_VARS_END
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        initData(i, j, k, var) = tileDesc%level
                    end do
                end do
            end do
        end do
    end associate
end subroutine Simulation_initBlock

