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

#include "Milhoja.h"

!> @details
!! A module that encapsulates the cpu_tf_hydro task function, which should
!! be passed to the Orchestration unit for applying the task function to
!! tile wrappers.
!!
!! This code should be written entirely by a Milhoja CODE GENERATOR.
module cpu_tf_hydro_mod
    implicit none
    private

    public :: cpu_tf_hydro_Fortran
    public :: cpu_tf_hydro_Cpp2C

    !!!!!----- INTERFACES TO C-LINKAGE C++ FUNCTIONS FOR TIME ADVANCE UNIT
   interface
      !> C++ task function that TimeAdvance passes to Orchestration unit
      subroutine cpu_tf_hydro_Cpp2C(C_threadIndex, C_dataItemPtr) bind(c)
         use iso_c_binding,     ONLY : C_PTR
         use milhoja_types_mod, ONLY : MILHOJA_INT
         integer(MILHOJA_INT), intent(IN), value :: C_threadIndex
         type(C_PTR),          intent(IN), value :: C_dataItemPtr
      end subroutine cpu_tf_hydro_Cpp2C
   end interface

contains

    !> @brief CPU-based variant of the hydro advance task function
    !!
    !! @details
    !! This routine is *not* a task function that is called directly by a pthread
    !! executing C++ code.  Rather the task function is the C++ code behind the
    !! cpu_tf_hydro_Cpp2C interface.  Despite this, this routine is still part of
    !! the Task Function/Data Packet Fortran/C++ interoperability design and
    !! should satisy that design's requirements.
    subroutine cpu_tf_hydro_Fortran(external_hydro_op1_dt, &
                                    external_hydro_op1_eosMode, &
                                    tile_deltas, &
                                    tile_hi, &
                                    tile_interior, &
                                    tile_lo, &
                                    CC_1, &
                                    scratch_hydro_op1_auxC, &
                                    scratch_hydro_op1_flX, &
                                    scratch_hydro_op1_flY, &
                                    scratch_hydro_op1_flZ, &
                                    lbdd_CC_1, &
                                    lbdd_scratch_hydro_op1_auxC, &
#if MILHOJA_NDIM == 1
                                    lbdd_scratch_hydro_op1_flX)
#elif MILHOJA_NDIM == 2
                                    lbdd_scratch_hydro_op1_flX, &
                                    lbdd_scratch_hydro_op1_flY)
#elif MILHOJA_NDIM == 3
                                    lbdd_scratch_hydro_op1_flX, &
                                    lbdd_scratch_hydro_op1_flY, &
                                    lbdd_scratch_hydro_op1_flZ)
#endif
        use Hydro_advanceSolution_variants_mod, ONLY : Hydro_computeSoundSpeed_block_cpu, &
                                                       Hydro_computeFluxes_X_block_cpu, &
                                                       Hydro_computeFluxes_Y_block_cpu, &
                                                       Hydro_computeFluxes_Z_block_cpu, &
                                                       Hydro_updateSolution_block_cpu
        use Eos_interface,                      ONLY : Eos_wrapped

        real,    intent(IN)            :: external_hydro_op1_dt
        integer, intent(IN)            :: external_hydro_op1_eosMode
        real,    intent(IN)            :: tile_deltas(1:MILHOJA_MDIM)
        integer, intent(IN)            :: tile_hi(1:MILHOJA_MDIM)
        integer, intent(IN)            :: tile_interior(1:2, 1:MILHOJA_MDIM)
        integer, intent(IN)            :: tile_lo(1:MILHOJA_MDIM)
        integer, intent(IN)            :: lbdd_CC_1(1:4)
        real,    intent(INOUT), target :: CC_1(:, :, :, :)
        integer, intent(IN)            :: lbdd_scratch_hydro_op1_auxC(1:3)
        real,    intent(OUT)           :: scratch_hydro_op1_auxC(:, :, :)
        integer, intent(IN)            :: lbdd_scratch_hydro_op1_flX(1:4)
#if MILHOJA_NDIM >= 2
        integer, intent(IN)            :: lbdd_scratch_hydro_op1_flY(1:4)
#endif
#if MILHOJA_NDIM == 3
        integer, intent(IN)            :: lbdd_scratch_hydro_op1_flZ(1:4)
#endif
        real,    intent(OUT)           :: scratch_hydro_op1_flX(:, :, :, :)
        real,    intent(OUT)           :: scratch_hydro_op1_flY(:, :, :, :)
        real,    intent(OUT)           :: scratch_hydro_op1_flZ(:, :, :, :)

        real, pointer :: CC_1_ptr(:, :, :, :)

        NULLIFY(CC_1_ptr)

        CALL Hydro_computeSoundSpeed_block_cpu(tile_lo, tile_hi, &
                                               CC_1, &
                                               lbdd_CC_1, &
                                               scratch_hydro_op1_auxC, &
                                               lbdd_scratch_hydro_op1_auxC)
        CALL Hydro_computeFluxes_X_block_cpu(external_hydro_op1_dt, &
                                             tile_lo, tile_hi, &
                                             tile_deltas, &
                                             CC_1, &
                                             lbdd_CC_1, &
                                             scratch_hydro_op1_auxC, &
                                             lbdd_scratch_hydro_op1_auxC, &
                                             scratch_hydro_op1_flX, &
                                             lbdd_scratch_hydro_op1_flX)
#if MILHOJA_NDIM >= 2
        CALL Hydro_computeFluxes_Y_block_cpu(external_hydro_op1_dt, &
                                             tile_lo, tile_hi, &
                                             tile_deltas, &
                                             CC_1, &
                                             lbdd_CC_1, &
                                             scratch_hydro_op1_auxC, &
                                             lbdd_scratch_hydro_op1_auxC, &
                                             scratch_hydro_op1_flY, &
                                             lbdd_scratch_hydro_op1_flY)
#endif
#if MILHOJA_NDIM == 3
        CALL Hydro_computeFluxes_Z_block_cpu(external_hydro_op1_dt, &
                                             tile_lo, tile_hi, &
                                             tile_deltas, &
                                             CC_1, &
                                             lbdd_CC_1, &
                                             scratch_hydro_op1_auxC, &
                                             lbdd_scratch_hydro_op1_auxC, &
                                             scratch_hydro_op1_flZ, &
                                             lbdd_scratch_hydro_op1_flZ)
#endif
        CALL Hydro_updateSolution_block_cpu(tile_lo, tile_hi, &
                                            scratch_hydro_op1_flX, &
                                            scratch_hydro_op1_flY, &
                                            scratch_hydro_op1_flZ, &
                                            lbdd_scratch_hydro_op1_flX, &
                                            CC_1, &
                                            lbdd_CC_1)

        ! Pointer is required for Eos_wrapped to work.  It will likely assume
        ! that the pointer has the correct indices set, so we have to set
        ! indices here.  EoS is always different...
        CC_1_ptr(lbdd_CC_1(1):, lbdd_CC_1(2):, lbdd_CC_1(3):, lbdd_CC_1(4):) => CC_1
        CALL Eos_wrapped(external_hydro_op1_eosMode, tile_interior, CC_1_ptr)
        NULLIFY(CC_1_ptr)
    end subroutine cpu_tf_hydro_Fortran

end module cpu_tf_hydro_mod
