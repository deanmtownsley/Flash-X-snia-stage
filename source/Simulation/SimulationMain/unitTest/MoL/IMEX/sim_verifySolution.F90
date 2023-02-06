
!!****if* source/Simulation/SimulationMain/unitTest/MoL/IMEX/sim_verifySolution
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
!!  sim_verifySolution
!!
!!
!! SYNOPSIS
!!  call sim_verifySolution(real,    intent(in)  :: t,
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
!!    valid    : Is this a valid solution
!!    maxError : The maximum error present in the solution that was utilized
!!               to determine if the solution was valid
!!
!!***
subroutine sim_verifySolution(t, valid, maxError)
   use Simulation_data, only: sim_A, sim_mu, sim_sigma, sim_beta, sim_alpha

   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator, &
                             Grid_getCellCoords
   use Grid_iterator, only: Grid_iterator_t
   use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "constants.h"

   implicit none

   real, intent(in) :: t
   logical, intent(out) :: valid
   real, intent(out) :: maxError

   real, dimension(:, :, :, :), pointer :: vars

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   integer :: i, j, k
   integer, dimension(LOW:HIGH, MDIM) :: lim

   real, allocatable :: x(:)
   real :: u_actual, u_err

   real, parameter :: errorTolerance = 2.5d-2

   nullify (vars)

   maxError = 0d0
   valid = .false.

   call Grid_getTileIterator(itor, LEAF)

   TileLoop: do
      if (.not. itor%isValid()) exit TileLoop

      call itor%currentTile(tileDesc)

      call tileDesc%getDataPtr(vars, CENTER)

      lim = tileDesc%limits

      allocate (x(tileDesc%limits(LOW, IAXIS):tileDesc%limits(HIGH, IAXIS)))
      x = 0d0
      call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                              tileDesc%limits(LOW, :), tileDesc%limits(HIGH, :), x)

      do k = lim(LOW, KAXIS), lim(HIGH, KAXIS)
         do j = lim(LOW, JAXIS), lim(HIGH, JAXIS)
            do i = lim(LOW, IAXIS), lim(HIGH, IAXIS)
               u_actual = exp(-sim_beta*t)*sim_A*exp(-0.5d0*((x(i) - sim_alpha*t - sim_mu)/sim_sigma)**2)

               u_err = abs((vars(U_VAR, i, j, k) - u_actual))

               maxError = max(maxError, u_err)
            end do ! i
         end do ! j
      end do ! k

      call tileDesc%releaseDataPtr(vars, CENTER)

      deallocate (x)

      call itor%next()
   end do TileLoop

   call Grid_releaseTileIterator(itor)

   if (maxError .lt. errorTolerance) valid = .true.
end subroutine sim_verifySolution
