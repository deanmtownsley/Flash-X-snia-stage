#include "constants.h"
#include "Simulation.h"

!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! A Milhoja-specific implementation of this routine.  Please refer
!! to the documentation in this routine's stub for general interface
!! information.
!!
!! @todo Should the deltas be cached by Grid_init?
subroutine Grid_getDeltas(level, deltas)
    use milhoja_types_mod,   ONLY : MILHOJA_INT, &
                                    MILHOJA_REAL
    use milhoja_grid_mod,    ONLY : milhoja_grid_getDeltas

    use gr_milhojaInterface, ONLY : gr_checkMilhojaError

    implicit none

    integer, intent(IN)  :: level
    real,    intent(OUT) :: deltas(1:MDIM)

    integer(MILHOJA_INT) :: MH_level
    real(MILHOJA_REAL)   :: MH_deltas(1:MDIM)
    integer(MILHOJA_INT) :: MH_ierr

    integer :: i

    ! Preset to the value that an element should have if its index is above NDIM
    deltas(:) = 0.0

    MH_level = INT(level, kind=MILHOJA_INT)

    CALL milhoja_grid_getDeltas(MH_level, MH_deltas, MH_ierr)
    CALL gr_checkMilhojaError("Grid_getDeltas", MH_ierr)

    do i = 1, NDIM
        deltas(i) = REAL(MH_deltas(i))
    end do
end subroutine Grid_getDeltas

