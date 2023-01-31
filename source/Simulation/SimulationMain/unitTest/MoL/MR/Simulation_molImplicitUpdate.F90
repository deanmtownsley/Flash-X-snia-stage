!!****if* source/Simulation/SimulationMain/unitTest/MoL/MR/Simulation_molImplicitUpdate
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
!!  NAME
!!
!!      Simulation_molImplicitUpdate
!!
!!  SYNOPSIS
!!
!!      call Simulation_molImplicitUpdate(real, intent(in) :: t
!!                                        real, intent(in) :: dt)
!!
!!  DESCRIPTION
!!
!!      Implicitly update evolved variables from t to t+dt
!!
!!
!!  ARGUMENTS
!!
!!      t  : Current time
!!      dt : Size of the time step to take
!!
!!***
subroutine Simulation_molImplicitUpdate(t, dt)
   use Simulation_data, only: sim_alpha, sim_beta, sim_epsilon, sim_lambdaF, sim_lambdaS

   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
   use Grid_iterator, only: Grid_iterator_t
   use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "constants.h"

   implicit none

   real, intent(in) :: t, dt

   integer, dimension(LOW:HIGH, MDIM) :: lim
   integer :: i, j, k, n

   real, dimension(:, :, :, :), pointer :: vars

   type(Grid_tile_t) :: tileDesc
   type(Grid_iterator_t) :: itor

   real :: u, v, A, B, cost, cosbt, f, dfdv, dv, v0

   integer, parameter :: max_iterations = 100
   real, parameter :: tol = 1d-10

   A = -0.5d0*sim_alpha*sim_epsilon*(sim_lambdaF - sim_lambdaS)
   B = 0.5d0*sim_lambdaS

   cost = cos(t)
   cosbt = cos(sim_beta*t)

   nullify (vars)

   call Grid_getTileIterator(itor, LEAF)

   do ! not using while since this is technically deprecated
      if (.not. itor%isValid()) exit

      call itor%currentTile(tileDesc)

      lim = tileDesc%limits

      ! This is equivalent to `call MoL_getDataPtr(tileDesc, vars, MOL_EVOLVED)`
      ! that is utilized in the RHS procedures of this simulation
      call tileDesc%getDataPtr(vars, CENTER)

      do k = lim(LOW, KAXIS), lim(HIGH, KAXIS)
         do j = lim(LOW, JAXIS), lim(HIGH, JAXIS)
            do i = lim(LOW, IAXIS), lim(HIGH, IAXIS)
               u = vars(U_VAR, i, j, k)
               v = vars(V_VAR, i, j, k)

               v0 = v

               ! Newton solve
               NewtonIteration: do n = 1, max_iterations
                  f = v - v0 - dt*(A*(u**2 - cosbt - 3d0)/u + B*(v**2 - cost - 2d0)/v)

                  ! if (f .eq. 0d0) exit NewtonIteration

                  dfdv = 1d0 - dt*(B*(v**2 + cost + 2d0)/v**2)
                  dv = f/dfdv

                  v = v - dv

                  if (abs(dv) .lt. tol) exit NewtonIteration
               end do NewtonIteration

               vars(V_VAR, i, j, k) = v
            end do ! i
         end do ! j
      end do ! k

      call tileDesc%releaseDataPtr(vars, CENTER)

      call itor%next()
   end do ! itor

   call Grid_releaseTileIterator(itor)

   call Simulation_molPostUpdate(t)

end subroutine Simulation_molImplicitUpdate
