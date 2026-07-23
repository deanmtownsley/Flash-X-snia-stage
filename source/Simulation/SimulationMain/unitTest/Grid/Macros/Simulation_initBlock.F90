!!****if* source/Simulation/SimulationMain/unitTest/Grid/Macros/Simulation_initBlock
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  call Simulation_initBlock(real,pointer :: initData(:,:,:,:),
!!                            integer(IN)  :: tileDesc   )
!!
!!
!!
!! DESCRIPTION
!!  This routine applies initial conditions of a specific simulation
!!  to the specified block.
!!
!!
!! ARGUMENTS
!!
!!  initData  -  pointer to solution data
!!  tileDesc  -  describes the block to initialize
!!
!!***

#include "constants.h"
#include "Simulation.h"

!!REORDER(4): initData
subroutine Simulation_initBlock(initData, tileDesc)
   use Grid_tile, ONLY : Grid_tile_t
   use Grid_interface, ONLY: Grid_getCellCoords

   implicit none

   real,              intent(IN), pointer :: initData(:, :, :, :)
   type(Grid_tile_t), intent(IN)          :: tileDesc

   real, allocatable :: xCenter(:)
   integer :: i, j, k

   associate(lo => tileDesc%limits(LOW,  :), &
         hi => tileDesc%limits(HIGH, :))
      allocate(xCenter(lo(IAXIS):hi(IAXIS)))
      call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xCenter)
      do k = lo(KAXIS), hi(KAXIS)
         do j = lo(JAXIS), hi(JAXIS)
            do i = lo(IAXIS), hi(IAXIS)
               initData(DENS_VAR, i, j, k) = 1.1*EXP(-xCenter(i)**2)
            end do
         end do
      end do
      deallocate(xCenter)
   end associate
end subroutine Simulation_initBlock

