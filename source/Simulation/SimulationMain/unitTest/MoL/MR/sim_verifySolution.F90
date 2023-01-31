subroutine sim_verifySolution(t, valid, maxError)
   use Simulation_data, only: sim_beta

   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
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

   real :: u_actual, v_actual, u_err, v_err

   real, parameter :: errorTolerance = 1d-6

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

               print *, i, j, k, u_err, v_err

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
end subroutine sim_verifySolution
