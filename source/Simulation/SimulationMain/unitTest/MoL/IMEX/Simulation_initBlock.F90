!!****if* source/Simulation/SimulationMain/unitTest/MoL/IMEX/Simulation_initBlock
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
subroutine Simulation_initBlock(vars, tileDesc)
   use Simulation_data, only: sim_A, sim_mu, sim_sigma

   use Grid_tile, only: Grid_tile_t
   use Grid_interface, only: Grid_getCellCoords

#include "Simulation.h"
#include "constants.h"

   implicit none

   real, dimension(:, :, :, :), pointer :: vars
   type(Grid_tile_t) :: tileDesc

   integer :: i, j, k
   real, allocatable :: x(:)

   allocate (x(tileDesc%limits(LOW, IAXIS):tileDesc%limits(HIGH, IAXIS)))
   x = 0d0
   call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                           tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :), x)

   do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
      do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
         do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
            vars(U_VAR, i, j, k) = sim_A*exp(-0.5d0*((x(i) - sim_mu)/sim_sigma)**2)
         end do ! i
      end do ! j
   end do ! k

   deallocate (x)

end subroutine Simulation_initBlock
