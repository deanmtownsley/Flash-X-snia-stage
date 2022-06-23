subroutine Simulation_initBlock(solnData, tileDesc)
    use Simulation_data

    use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "constants.h"

    implicit none

    real, dimension(:,:,:,:), pointer :: solnData
    type(Grid_tile_t) :: tileDesc

    real :: del(MDIM), box(LOW:HIGH,MDIM)
    real :: x
    integer :: i, j, k

    call tileDesc%deltas(del)
    call tileDesc%boundBox(box)

    do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
        do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
            do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                x = box(LOW,IAXIS) + (i - 0.5d0)*del(IAXIS)
                solnData((/PHI0_VAR,PHI1_VAR/),i,j,k) = &
                    (/1d0,2d0/)*sim_height*exp(-0.5d0*((x - sim_center)/sim_width)**2)
            end do ! i
        end do ! j
    end do ! k

end subroutine Simulation_initBlock
