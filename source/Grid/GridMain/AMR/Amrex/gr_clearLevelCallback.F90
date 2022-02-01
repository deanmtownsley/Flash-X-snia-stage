#ifdef DEBUG_ALL
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
#define DEBUG_GRID
#endif

#include "Simulation.h"

subroutine gr_clearLevelCallback(lev) bind(c)
    use amrex_amr_module,          ONLY : amrex_multifab_destroy
    use gr_fluxregister_mod,       ONLY : gr_fluxregisterDestroy

    use Grid_data,                 ONLY : gr_doFluxCorrection, &
                                          gr_amrexDidRefinement
    use gr_physicalMultifabs,      ONLY : unk, &
                                          gr_scratchCtr, &
                                          facevarx, facevary, facevarz, &
                                          fluxes, &
                                          flux_registers

    implicit none

    integer, intent(in), value :: lev

    integer :: dir

    ! Communicate to Grid_updateRefinement that this level might have been
    ! completely derefined
    gr_amrexDidRefinement = .TRUE.

    ! Multifab arrays use 0-based index set like AMReX
    call amrex_multifab_destroy(unk     (lev))
#if NFACE_VARS > 0
    call amrex_multifab_destroy(facevarx(lev))
#if NDIM >= 2
    call amrex_multifab_destroy(facevary(lev))
#endif
#if NDIM == 3
    call amrex_multifab_destroy(facevarz(lev))
#endif
#endif

    if (allocated(gr_scratchCtr))  call amrex_multifab_destroy(gr_scratchCtr(lev))

#if NFLUXES > 0
# ifdef USE_LEVELWIDE_FLUXES
    do dir = 1, SIZE(fluxes, 2)
        call amrex_multifab_destroy(fluxes(lev, dir))
    end do
# endif

    if ((lev > 0) .AND. (gr_doFluxCorrection)) then
        call gr_fluxregisterDestroy(flux_registers(lev))
    end if
#endif

#ifdef DEBUG_GRID
    write(*,'(A,A,I2)') "[gr_clearLevelCallback]", &
                        "              Cleared level", lev + 1
#endif

end subroutine gr_clearLevelCallback

