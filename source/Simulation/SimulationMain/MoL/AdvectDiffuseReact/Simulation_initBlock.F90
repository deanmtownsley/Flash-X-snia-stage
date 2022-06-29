subroutine Simulation_initBlock(vars, tileDesc)
    use Simulation_data
    use Driver_interface, only: Driver_getSimTime
    use Grid_tile, only: Grid_tile_t
    use Grid_interface, only: Grid_getCellCoords

#include "Simulation.h"
#include "constants.h"

    implicit none

    real, dimension(:,:,:,:), pointer :: vars
    type(Grid_tile_t) :: tileDesc

    real :: del(MDIM), box(LOW:HIGH,MDIM)
    integer :: i, j, k
    real :: t

    real, allocatable :: x(:)

    call Driver_getSimTime(t)

    call tileDesc%deltas(del)
    call tileDesc%boundBox(box)

    allocate(x(tileDesc%limits(LOW,IAXIS):tileDesc%limits(HIGH,IAXIS)))
    x = 0d0
    call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                            tileDesc%limits(LOW,:), tileDesc%limits(HIGH,:), &
                            x)

    do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
        do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
            do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                vars(U_VAR,i,j,k) = sim_a       + 0.1d0*sin(PI*x(i))
                vars(V_VAR,i,j,k) = sim_b/sim_a + 0.1d0*sin(PI*x(i))
                vars(W_VAR,i,j,k) = sim_b       + 0.1d0*sin(PI*x(i))
            end do ! i
        end do ! j
    end do ! k

    deallocate(x)

end subroutine Simulation_initBlock
