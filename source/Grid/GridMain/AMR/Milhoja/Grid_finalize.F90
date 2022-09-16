#include "constants.h"

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
subroutine Grid_finalize()
    use milhoja_types_mod,   ONLY : MILHOJA_INT
    use milhoja_grid_mod,    ONLY : milhoja_grid_finalize

    use Grid_data,           ONLY : gr_globalMe
    use gr_milhojaInterface, ONLY : gr_checkMilhojaError

    implicit none

    integer(MILHOJA_INT) :: MH_ierr

    if (gr_globalMe == MASTER_PE) then
        write(*,'(A)') "[Grid] Finalizing ..."
    end if

    CALL milhoja_grid_finalize(MH_ierr)
    CALL gr_checkMilhojaError("Grid_finalize", MH_ierr)

    if (gr_globalMe == MASTER_PE) then
        write(*,'(A)') "[Grid] Finalized"
    end if
end subroutine Grid_finalize

