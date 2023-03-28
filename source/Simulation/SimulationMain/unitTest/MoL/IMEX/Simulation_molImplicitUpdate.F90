!!****if* source/Simulation/SimulationMain/unitTest/MoL/IMEX/Simulation_molImplicitUpdate
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

!!REORDER(4): vars

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

   real :: u, v, Au, Bu, Av, Bv, cost, cosbt
   real :: f(2), dfdv(2, 2), invdfdv(2, 2), detdfdv, du, dv, u0, v0

   integer, parameter :: max_iterations = 100
   real, parameter :: tol = 1e-10

   Au = 0.5*sim_lambdaF
   Bu = 0.5*(1.0 - sim_epsilon)*(sim_lambdaF - sim_lambdaS)/sim_alpha

   Av = -0.5*sim_alpha*sim_epsilon*(sim_lambdaF - sim_lambdaS)
   Bv = 0.5*sim_lambdaS

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

               u0 = u
               v0 = v

               ! Newton solve
               NewtonIteration: do n = 1, max_iterations
                  f(1) = u - u0 - dt*((Au*(u**2 - cosbt - 3.0) - 0.5*sim_beta*sin(sim_beta*t))/u &
                                      + Bu*(v**2 - cost - 2.0)/v)
                  f(2) = v - v0 - dt*(Av*(u**2 - cosbt - 3.0)/u + Bv*(v**2 - cost - 2.0)/v)

                  ! if (f .eq. 0.0) exit NewtonIteration

                  dfdv(1, 1) = 1.0 - dt*(2.0*Au*(u**2 + cosbt + 3.0) + 0.5*sim_beta*sin(sim_beta*t))/(2.0*u**2)
                  dfdv(1, 2) = -dt*Bu*(v**2 + cost + 2.0)/v**2
                  dfdv(2, 1) = -dt*(Av*(u**2 + cosbt + 3.0)/u**2)
                  dfdv(2, 2) = 1.0 - dt*(Bv*(v**2 + cost + 2.0)/v**2)

                  detdfdv = dfdv(1, 1)*dfdv(2, 2) - dfdv(1, 2)*dfdv(2, 1)

                  invdfdv(1, 1) = dfdv(2, 2)/detdfdv
                  invdfdv(1, 2) = -dfdv(2, 1)/detdfdv
                  invdfdv(2, 1) = -dfdv(1, 2)/detdfdv
                  invdfdv(2, 2) = dfdv(1, 1)/detdfdv

                  du = invdfdv(1, 1)*f(1) + invdfdv(1, 2)*f(2)
                  dv = invdfdv(2, 1)*f(1) + invdfdv(2, 2)*f(2)

                  u = u - du
                  v = v - dv

                  if (sqrt(du**2 + dv**2) .lt. tol) exit NewtonIteration
               end do NewtonIteration

               vars(U_VAR, i, j, k) = u
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
