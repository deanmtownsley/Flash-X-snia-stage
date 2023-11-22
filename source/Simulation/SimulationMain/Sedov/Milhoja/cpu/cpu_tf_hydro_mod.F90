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
!! @file

#include "constants.h"
#include "Simulation.h"

! The use of this module implies that the runtime is being used.
! Therefore, the Driver unit can assume that precision checks were
! performed by Orchestration_init.  Hence, type casts need not be checked.
#ifndef FLASHX_ORCHESTRATION_MILHOJA
#error "Task function bundles imply the use of the Milhoja runtime"
#endif

!> @details
!! A module that encapsulates the cpu_tf_hydro task function, which should
!! passed to the Orchestration unit for applying the task function to tiles.
!!
!! This code should be written entirely by a Milhoja CODE GENERATOR as part of
!! the full set of bundles needed to execute the timestep recipe.
!!
!! @todo This is presently Flash-X specific due to the use of Grid_tile_t and
!! the inclusion of Flash-X header files, which means that Milhoja cannot
!! generate this yet.  Update this scheme so that it matches what Wesley has
!! done with the data packet.
module cpu_tf_hydro_mod
    implicit none
    private

    public :: cpu_tf_hydro

contains

    !> @brief CPU-based variant of the hydro advance task function
    !!
    !! @details
    !! Since this routine is called directly by a pthread running C++ code, it
    !! is part of the Task Function/Data Packet Fortran/C++ interoperability
    !! design and should satisy that design's requirements.  This justifies the
    !! fact that the interface exposes iso_c_binding types.
    subroutine cpu_tf_hydro(C_threadId, C_dataItemPtr) bind(c)
        use iso_c_binding, ONLY : C_PTR, &
                                  C_NULL_PTR, &
                                  C_F_POINTER

        use milhoja_types_mod,                  ONLY : MILHOJA_INT, &
                                                       MILHOJA_REAL
        use milhoja_tile_mod,                   ONLY : milhoja_tile_from_wrapper_C
        use cpu_tf_hydro_C2F_mod,               ONLY : get_dt_wrapper_C, &
                                                       get_scratch_auxC_wrapper_C

        use Grid_tile,                          ONLY : Grid_tile_t, &
                                                       Grid_tile_fromMilhojaTilePtr
        use gr_milhojaInterface,                ONLY : gr_checkMilhojaError
        use Eos_interface,                      ONLY : Eos_wrapped
        use Hydro_data,                         ONLY : hy_eosModeAfter
        use Hydro_advanceSolution_variants_mod, ONLY : Hydro_computeSoundSpeed_block_cpu, &
                                                       Hydro_computeFluxes_X_block_cpu, &
                                                       Hydro_computeFluxes_Y_block_cpu, &
                                                       Hydro_computeFluxes_Z_block_cpu, &
                                                       Hydro_updateSolution_block_cpu

        integer(MILHOJA_INT), intent(IN), value :: C_threadId
        type(C_PTR),          intent(IN), value :: C_dataItemPtr

        integer(MILHOJA_INT)         :: MH_ierr
        integer                      :: F_threadId
        type(C_PTR)                  :: C_tilePtr
        type(Grid_tile_t)            :: F_tile
        type(C_PTR)                  :: C_auxC
        real,                pointer :: F_auxC(:, :, :)
        integer                      :: lbdd_auxC(1:3)
        real(MILHOJA_REAL)           :: MH_dt
        real                         :: F_dt

        real, pointer                    :: U(:, :, :, :)
        real, pointer                    :: flX(:, :, :, :)
        real, pointer                    :: flY(:, :, :, :)
        real, pointer                    :: flZ(:, :, :, :)
        real,                     target :: empty4(0, 0, 0, 0)
        real                             :: deltas(1:MDIM)

        C_tilePtr = C_NULL_PTR
        C_auxC = C_NULL_PTR
        NULLIFY(F_auxC)
        NULLIFY(U)
        NULLIFY(flX)
        NULLIFY(flY)
        NULLIFY(flZ)

        !!!!!----- CONVERT TileWrapper pointer to Tile pointer
        MH_ierr = get_dt_wrapper_C(C_dataItemPtr, MH_dt)
        CALL gr_checkMilhojaError("cpu_tf_hydro_mod", MH_ierr)
        F_dt = REAL(MH_dt)

        MH_ierr = get_scratch_auxC_wrapper_C(C_dataItemPtr, &
                                             C_threadId, &
                                             C_auxC)
        CALL gr_checkMilhojaError("cpu_tf_hydro_mod", MH_ierr)
        CALL C_F_POINTER(C_auxC, F_auxC, shape=[18, 18, 18])
        C_auxC = C_NULL_PTR

        MH_ierr = milhoja_tile_from_wrapper_C(C_dataItemPtr, C_tilePtr)
        CALL gr_checkMilhojaError("cpu_tf_hydro_mod", MH_ierr)

        !!!!!----- CONVERT C-ARGUMENTS TO FORTRAN VARIABLES
        F_threadId = INT(C_threadId)
        CALL Grid_tile_fromMilhojaTilePtr(C_tilePtr, F_tile)

        !!!!!----- PULL APART DATA ITEM
        ! Map contents of dataItem onto argument list of static Fortran routines
        CALL F_tile%deltas(deltas)
        CALL F_tile%getDataPtr(U,   CENTER)
        CALL F_tile%getDataPtr(flX, FLUXX)
        CALL F_tile%getDataPtr(flY, FLUXY)
        CALL F_tile%getDataPtr(flZ, FLUXZ)
        if (.NOT. ASSOCIATED(flY)) then
            flY => empty4
        end if
        if (.NOT. ASSOCIATED(flZ)) then
            flZ => empty4
        end if

        !!!!!----- ADVANCE SOLUTION IN INTERIOR
        associate(lo => F_tile%limits(LOW,  :), &
                  hi => F_tile%limits(HIGH, :))
            lbdd_auxC = [lo(IAXIS) - K1D, lo(JAXIS) - K2D, lo(KAXIS) - K3D]

            CALL Hydro_computeSoundSpeed_block_cpu(lo, hi, &
                                                   U,    lbound(U), &
                                                   F_auxC, lbdd_auxC)
            CALL Hydro_computeFluxes_X_block_cpu(F_dt, lo, hi, deltas, &
                                                 U,    lbound(U), &
                                                 F_auxC, lbdd_auxC, &
                                                 flX,  lbound(flX))
#if NDIM >= 2
            CALL Hydro_computeFluxes_Y_block_cpu(F_dt, lo, hi, deltas, &
                                                 U,    lbound(U), &
                                                 F_auxC, lbdd_auxC, &
                                                 flY,  lbound(flY))
#endif
#if NDIM == 3
            CALL Hydro_computeFluxes_Z_block_cpu(F_dt, lo, hi, deltas, &
                                                 U,    lbound(U), &
                                                 F_auxC, lbdd_auxC, &
                                                 flZ,  lbound(flZ))
#endif
            NULLIFY(F_auxC)
            CALL Hydro_updateSolution_block_cpu(lo, hi, &
                                                flX, flY, flZ, lbound(flX), &
                                                U, lbound(U))
            CALL Eos_wrapped(hy_eosModeAfter, F_tile%limits, U)
        end associate

        CALL F_tile%releaseDataPtr(flZ, FLUXZ)
        CALL F_tile%releaseDataPtr(flY, FLUXY)
        CALL F_tile%releaseDataPtr(flX, FLUXX)
        CALL F_tile%releaseDataPtr(U,   CENTER)
    end subroutine cpu_tf_hydro

end module cpu_tf_hydro_mod

