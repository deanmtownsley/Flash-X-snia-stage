!!****if* source/Grid/GridMain/AMR/Amrex/Grid_communicateFluxes
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
!!  Grid_communicateFluxes
!!
!! SYNOPSIS
!!  call Grid_communicateFluxes(integer(IN) :: axis,
!!                              integer(IN) :: coarse_level)
!!
!! DESCRIPTION
!!
!!  This call launches the communication phase for flux correction.
!!
!!  The communication happens entirely within SPFS, and does not
!!  involve any data structures that should be directly visible to
!!  client code. Instead, client code uses subroutines like
!!  Grid_putFluxData_fluxbuf to transfer flux data to SPFS before a
!!  communication phase, and Grid_correctFluxData to get data back
!!  from SPFS after the communication.
!!
!!  The coarse_level argument selects a pair of refinement levels
!!  whose blocks may participate in the flux correction.  The
!!  coarse_level argument identifies the refinement level on the
!!  coarse size of the fine/coarse boundaries for which data is
!!  exchanged.
!!
!!  It is assumed that before calling this routine, the code has
!!  already stored the correct fluxes for all relevant blocks at the
!!  finer level into the semi-permanent flux storage, using
!!  Grid_putFluxData or a variant thereof.
!!
!!  This routine will comunicate the data necessary, including
!!  communication between MPI ranks, so that the flux data on the
!!  coarse side of fine/coarse boundaries can be overwritten by
!!  presumably more accurate values from computations on the fine
!!  side.
!!
!! ARGUMENTS 
!!  axis - the only acceptable value for AMReX is ALLDIR.
!!  coarse_level - the 1-based level index of the coarse blocks
!!
!! NOTES
!!
!!   SPFS means semi-permanent flux storage. When using a Grid
!!   implementation based on AMReX, SPFS is implemented by an AMReX
!!   flux register class, such as FlashFluxRegister.
!!
!!   The special value UNSPEC_LEVEL for the coarse_value argument,
!!   which should be taken to mean "communicate between all levels",
!!   is not supported in this implementation for AMReX.
!!   That is, a specific (coarse) level is always required.
!!
!! SEE ALSO
!!
!!  Grid_putFluxData
!!  Grid_putFluxData_fluxbuf
!!  Grid_correctFluxData
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "Simulation.h"
#include "constants.h"

subroutine Grid_communicateFluxes(axis, coarse_level)
    use gr_physicalMultifabs, ONLY : flux_registers
    use amrex_amrcore_module, ONLY : amrex_get_finest_level

    implicit none

    integer, intent(IN)                   :: axis
    integer, intent(IN)                   :: coarse_level

    integer :: fine
    integer :: coarse

    if (axis /= ALLDIR) then
        call Driver_abort("[Grid_conserveFluxes] AMReX currently requires axis==ALLDIR")
    end if

#ifndef USE_AMREX_FLASHFLUXREGISTER
    call Driver_abort("Grid_communicateFluxes.F90 requires amrex_flash_fluxregister,&
       & make sure USE_AMREX_FLASHFLUXREGISTER is defined!")
#else

    if (coarse_level == UNSPEC_LEVEL) then
       call Driver_abort("This Grid_communicateFluxes.F90 does not support UNSPEC_LEVEL!")
    end if

    ! FLASH uses 1-based level index / AMReX uses 0-based index
    coarse = coarse_level - 1
    fine   = coarse_level

    ! No need to communicate on the finest level in existence or any
    ! level index corresponding to a finer mesh
    if (coarse >= amrex_get_finest_level())     RETURN

        ! The AMReX flux registers internally are storing and communicating
        ! fluxes and *not* flux densities.

    call flux_registers(fine)%communicate()

#endif
end subroutine Grid_communicateFluxes

