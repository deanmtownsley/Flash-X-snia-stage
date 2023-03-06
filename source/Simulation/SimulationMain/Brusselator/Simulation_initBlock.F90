!!****if* source/Simulation/SimulationMain/Brusselator/Simulation_initBlock
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
!!  call Simulation_initBlock(real,pointer :: solnData(:,:,:,:),
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
!!  solnData  -        pointer to solution data
!!  blockDesc -        describes the block to initialize
!!
!!***

!!REORDER(4): vars

subroutine Simulation_initBlock(vars, tileDesc)
   use Simulation_data, only: sim_b, sim_a

   use Driver_interface, only: Driver_getSimTime
   use Grid_tile, only: Grid_tile_t
   use Grid_interface, only: Grid_getCellCoords

#include "Simulation.h"
#include "constants.h"

   implicit none

   real, dimension(:, :, :, :), pointer :: vars
   type(Grid_tile_t) :: tileDesc

   real :: del(MDIM), box(LOW:HIGH, MDIM)
   integer :: i, j, k
   real :: t

   real, allocatable :: x(:)

   call Driver_getSimTime(t)

   call tileDesc%deltas(del)
   call tileDesc%boundBox(box)

   allocate (x(tileDesc%limits(LOW, IAXIS):tileDesc%limits(HIGH, IAXIS)))
   x = 0.0
   call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                           tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :), x)

   do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
      do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
         do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
            vars(U_VAR, i, j, k) = sim_a + 0.1*sin(PI*x(i))
            vars(V_VAR, i, j, k) = sim_b/sim_a + 0.1*sin(PI*x(i))
            vars(W_VAR, i, j, k) = sim_b + 0.1*sin(PI*x(i))
         end do ! i
      end do ! j
   end do ! k

   deallocate (x)

end subroutine Simulation_initBlock
