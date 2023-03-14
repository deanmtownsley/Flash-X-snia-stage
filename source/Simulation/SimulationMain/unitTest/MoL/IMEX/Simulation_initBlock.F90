!!****if* source/Simulation/SimulationMain/unitTest/MoL/IMEX/Simulation_initBlock
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
!!  call Simulation_initBlock(real,pointer :: vars(:,:,:,:),
!!                            integer(IN)  :: blockDesc  )
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
!!  vars      -        pointer to solution data
!!  blockDesc -        describes the block to initialize
!!
!!***

!!REORDER(4): vars

subroutine Simulation_initBlock(vars, tileDesc)
   use Simulation_data, only: sim_beta

   use Driver_interface, only: Driver_getSimTime
   use Grid_tile, only: Grid_tile_t
   use Grid_interface, only: Grid_getCellCoords

#include "Simulation.h"
#include "constants.h"

   implicit none

   real, dimension(:, :, :, :), pointer :: vars
   type(Grid_tile_t) :: tileDesc

   integer :: i, j, k
   real :: t

   call Driver_getSimTime(t)

   do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
      do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
         do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
            vars(U_VAR, i, j, k) = 2.0
            vars(V_VAR, i, j, k) = sqrt(3.0)
         end do ! i
      end do ! j
   end do ! k

end subroutine Simulation_initBlock
