!!****if* source/Simulation/SimulationMain/Brusselator/Simulation_initBlock
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
subroutine Simulation_initBlock(vars, tileDesc)
   !!use Simulation_data, only: sim_b, sim_a

   use Driver_interface, only: Driver_getSimTime
   use Grid_tile, only: Grid_tile_t
   use Grid_interface, only: Grid_getCellCoords

#include "Z4c.h"
#include "constants.h"

   implicit none

   real, dimension(:, :, :, :), pointer :: vars
   type(Grid_tile_t) :: tileDesc

   real :: del(MDIM), box(LOW:HIGH, MDIM)
   integer :: i, j, k
   real :: t

   real :: A, d, b, dtb, dxb, dyb, dzb, chi, Ksclr, &
           gamma_LL_00, gamma_LL_11, gamma_LL_22, K_LL_00, K_LL_11, K_LL_22, x_U_0, x_U_1, x_U_2

   real, allocatable :: x(:)
   real, allocatable :: y(:)
   real, allocatable :: z(:)

   call Driver_getSimTime(t)

   call tileDesc%deltas(del)
   call tileDesc%boundBox(box)

   allocate (x(tileDesc%limits(LOW, IAXIS):tileDesc%limits(HIGH, IAXIS)))
   x = 0d0
   call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                           tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :), x)

   allocate (y(tileDesc%limits(LOW, JAXIS):tileDesc%limits(HIGH, JAXIS)))
   y = 0d0
   call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                           tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :), y)

   allocate (z(tileDesc%limits(LOW, KAXIS):tileDesc%limits(HIGH, KAXIS)))
   z = 0d0
   call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                           tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :), z)

   do k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
      do j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
         do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)

            x_U_0 = x(i)
            x_U_1 = y(j)
            x_U_2 = z(k)

            A = 0.05d0
            d = 1d0
            b = A*sin(2d0*PI*x_U_0/d)
            dtb = -(2d0*PI*A/d)*cos(2*PI*x_U_0/d)
            dxb = (2d0*PI*A/d)*cos(2*PI*x_U_0/d)
            chi = (1d0/(1d0 - b))**(1.0d0/3.0d0)

            Ksclr = (dtb/2d0)*(1d0/(1d0 - b))**1.5d0

            gamma_LL_00 = 1d0 - b
            gamma_LL_11 = 1d0
            gamma_LL_22 = 1d0
            K_LL_00 = (dtb/2d0)*sqrt(1d0/(1d0 - b))
            K_LL_11 = 0
            K_LL_22 = 0

            vars(CHI_VAR, i, j, k) = chi

            vars(GAMTILDE_LL_00_VAR, i, j, k) = chi*(1d0 - b)
            vars(GAMTILDE_LL_01_VAR, i, j, k) = 0d0
            vars(GAMTILDE_LL_02_VAR, i, j, k) = 0d0
            vars(GAMTILDE_LL_11_VAR, i, j, k) = chi
            vars(GAMTILDE_LL_12_VAR, i, j, k) = 0d0
            vars(GAMTILDE_LL_22_VAR, i, j, k) = chi

            vars(KHAT_VAR, i, j, k) = Ksclr

            vars(ATILDE_LL_00_VAR, i, j, k) = chi*(K_LL_00 - (1d0/3d0)*gamma_LL_00*Ksclr)
            vars(ATILDE_LL_01_VAR, i, j, k) = 0d0
            vars(ATILDE_LL_02_VAR, i, j, k) = 0d0
            vars(ATILDE_LL_11_VAR, i, j, k) = chi*(K_LL_11 - (1d0/3d0)*gamma_LL_11*Ksclr)
            vars(ATILDE_LL_12_VAR, i, j, k) = 0d0
            vars(ATILDE_LL_22_VAR, i, j, k) = chi*(K_LL_22 - (1d0/3d0)*gamma_LL_22*Ksclr)

            vars(GAMTILDE_U_0_VAR, i, j, k) = -(2d0*dxb/3d0)*(1d0/(1d0 - b))**(5d0/3d0)
            vars(GAMTILDE_U_1_VAR, i, j, k) = 0d0
            vars(GAMTILDE_U_2_VAR, i, j, k) = 0d0

            vars(ALPHA_VAR, i, j, k) = sqrt(1d0 - b)

            vars(BETA_U_0_VAR, i, j, k) = 0d0
            vars(BETA_U_1_VAR, i, j, k) = 0d0
            vars(BETA_U_2_VAR, i, j, k) = 0d0
         end do ! i
      end do ! j
   end do ! k

   deallocate (x)
   deallocate (y)
   deallocate (z)

end subroutine Simulation_initBlock
