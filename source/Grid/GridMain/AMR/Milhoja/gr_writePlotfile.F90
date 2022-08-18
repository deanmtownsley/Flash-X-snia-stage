#include "Milhoja_interface_error_codes.h"

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
!! A Milhoja-specific routine that requests that the Milhoja grid backend dump the
!! contents of its Grid data structures to a plotfile.  The data is written to a
!! file with name
!!                          milhoja_plt_<step> 
!!
!! It is intended that this only be used for development, testing, and debugging.
!! Therefore, this routine has not been added to the public interface and it is
!! reasonable for non-Grid code (e.g, Driver) to call on a temporary basis this
!! Grid-private routine in violation of the Flash-X design rules.
!!
!! @param step   The number of the current time step
subroutine gr_writePlotfile(step)
    use mpi

    use milhoja_types_mod,   ONLY : MILHOJA_INT
    use milhoja_grid_mod,    ONLY : milhoja_grid_writePlotfile

    use Grid_data,           ONLY : gr_globalComm
    use gr_milhojaInterface, ONLY : gr_checkMilhojaError

    implicit none

    integer, intent(IN) :: step

    integer(MILHOJA_INT) :: MH_step
    integer(MILHOJA_INT) :: MH_ierr

    integer :: mpierr

    CALL MPI_Barrier(gr_globalComm, mpierr)

    MH_step = INT(step, kind=MILHOJA_INT)
    CALL milhoja_grid_writePlotfile(MH_step, MH_ierr)
    CALL gr_checkMilhojaError("gr_writePlotfile", MH_ierr)
end subroutine gr_writePlotfile

