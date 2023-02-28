
!!****if* source/Simulation/SimulationMain/unitTest/MoL/IMEX/sim_molVerifySolution
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
!!  sim_molVerifySolution
!!
!!
!! SYNOPSIS
!!  call sim_molVerifySolution(real,    intent(in)  :: t,
!!                          real,    intent(in)  :: dt,
!!                          logical, intent(out) :: valid,
!!                          real,    intent(out) :: maxError)
!!
!!
!! DESCRIPTION
!!    Verify the solution at the specified time.  An implementation
!!    may utilize an analytic/baseline solution or some other
!!    error-estimator to determine if the current solution is valid
!!    for the purposes of passing a unit test.
!!
!! ARGUMENTS
!!
!!    t        : The current time that the solution is at
!!    dt       : The current timestep
!!    valid    : Is this a valid solution
!!    maxError : The maximum error present in the solution that was utilized
!!               to determine if the solution was valid
!!
!!***
subroutine sim_molVerifySolution(t, dt, valid, maxError)
   use Simulation_data, only: sim_beta

   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
   use Grid_iterator, only: Grid_iterator_t
   use Grid_tile, only: Grid_tile_t

   use MoL_interface, only: MoL_getOrder

#include "Simulation.h"
#include "constants.h"

   implicit none

   real, intent(in) :: t, dt
   logical, intent(out) :: valid
   real, intent(out) :: maxError

   real, dimension(:, :, :, :), pointer :: vars

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   integer :: i, j, k
   integer, dimension(LOW:HIGH, MDIM) :: lim

   real :: u_actual, v_actual, u_err, v_err

   real :: errorTolerance

   errorTolerance = 2d0*dt**(MoL_getOrder() - 1)

   nullify (vars)

   maxError = 0d0
   valid = .false.

   ! No guard-cell filling necessary - just a bunch of local equations to solve

   call Grid_getTileIterator(itor, LEAF)

   TileLoop: do
      if (.not. itor%isValid()) exit TileLoop

      call itor%currentTile(tileDesc)

      call tileDesc%getDataPtr(vars, CENTER)

      lim = tileDesc%limits

      do k = lim(LOW, KAXIS), lim(HIGH, KAXIS)
         do j = lim(LOW, JAXIS), lim(HIGH, JAXIS)
            do i = lim(LOW, IAXIS), lim(HIGH, IAXIS)
               u_actual = sqrt(3d0 + cos(sim_beta*t))
               v_actual = sqrt(2d0 + cos(t))

               u_err = abs(vars(U_VAR, i, j, k) - u_actual)
               v_err = abs(vars(V_VAR, i, j, k) - v_actual)

               maxError = max(maxError, u_err)
               maxError = max(maxError, v_err)
            end do ! i
         end do ! j
      end do ! k

      call tileDesc%releaseDataPtr(vars, CENTER)

      call itor%next()
   end do TileLoop

   call Grid_releaseTileIterator(itor)

   if (maxError .lt. errorTolerance) valid = .true.
end subroutine sim_molVerifySolution
