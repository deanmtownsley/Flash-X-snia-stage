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
!! Only partial functionality has been implemented so far.  This routine
!! aborts if calling code attempts to use non-implemented functionality.
!!
!! @todo Pass MH_T_INIT all the way through to the Grid backend.
!! @todo Do we use RPs to determine which initDomain routine to
!!       call (e.g., no runtime, CPU-only, GPU-only, etc.)
!! @todo Code up full implementation
subroutine Grid_initDomain(restart, particlesInitialized)
    use milhoja_types_mod,           ONLY : MILHOJA_INT, &
                                            MILHOJA_REAL
    use milhoja_grid_mod,            ONLY : milhoja_grid_initDomain

    use gr_milhojaInterface,         ONLY : gr_checkMilhojaError
    use gr_initDomain_mod,           ONLY : gr_initBlock_tile_cpu
    use Driver_interface,            ONLY : Driver_abort

    implicit none

    real(MILHOJA_REAL), parameter :: MH_T_INIT = 0.0

    logical, intent(IN)    :: restart
    logical, intent(INOUT) :: particlesInitialized

    integer(MILHOJA_INT) :: MH_ierr

    particlesInitialized = .FALSE.

    if (restart) then
        CALL Driver_abort("[Grid_initDomain] Restarts not yet implemented")
    end if

#ifdef USE_MILHOJA_RUNTIME
    CALL Driver_abort("[Grid_initDomain] Runtime not yet integrated")
#else
    CALL milhoja_grid_initDomain(gr_initBlock_tile_cpu, MH_ierr)
#endif
    CALL gr_checkMilhojaError("Grid_initDomain", MH_ierr)
end subroutine Grid_initDomain

