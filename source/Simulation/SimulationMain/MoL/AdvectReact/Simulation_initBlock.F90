subroutine Simulation_initBlock(solnData, tileDesc)
    use Simulation_data

    use Grid_tile, only: Grid_tile_t
    use Grid_interface, only: Grid_getCellCoords

#include "Simulation.h"
#include "constants.h"

    implicit none

    real, dimension(:,:,:,:), pointer :: solnData
    type(Grid_tile_t) :: tileDesc

    real :: del(MDIM), box(LOW:HIGH,MDIM)
    ! real :: x
    integer :: i, j, k

    real, allocatable :: x(:)

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
                ! x = box(LOW,IAXIS) + (i - 0.5d0)*del(IAXIS)
                solnData((/PHI0_VAR,PHI1_VAR/),i,j,k) = &
                    (/1d0,2d0/)*sim_height*exp(-0.5d0*((x(i) - sim_center)/sim_width)**2)

                solnData(DENS_VAR,i,j,k) = 1d0
            end do ! i
        end do ! j
    end do ! k

    deallocate(x)

end subroutine Simulation_initBlock
