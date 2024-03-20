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

subroutine cpu_tf_hydro_C2F(MH_external_hydro_op1_dt, &
                            MH_external_hydro_op1_eosMode, &
                            C_tile_deltas, &
                            C_tile_hi, &
                            C_tile_interior, &
                            C_tile_lo, &
                            C_CC_1, &
                            C_scratch_hydro_op1_auxC, &
                            C_scratch_hydro_op1_flX, &
                            C_scratch_hydro_op1_flY, &
                            C_scratch_hydro_op1_flZ, &
                            C_lbdd_CC_1, &
                            C_lbdd_scratch_hydro_op1_auxC, &
                            C_lbdd_scratch_hydro_op1_flX, &
                            C_lbdd_scratch_hydro_op1_flY) bind(c)
    use iso_c_binding, ONLY : C_PTR, &
                              C_F_POINTER

    use milhoja_types_mod, ONLY : MILHOJA_INT, &
                                  MILHOJA_REAL
    use cpu_tf_hydro_mod,  ONLY : cpu_tf_hydro_Fortran

    implicit none

    real(MILHOJA_REAL),   intent(IN), value :: MH_external_hydro_op1_dt
    integer(MILHOJA_INT), intent(IN), value :: MH_external_hydro_op1_eosMode
    type(C_PTR),          intent(IN), value :: C_tile_deltas
    type(C_PTR),          intent(IN), value :: C_tile_hi
    type(C_PTR),          intent(IN), value :: C_tile_interior
    type(C_PTR),          intent(IN), value :: C_tile_lo
    type(C_PTR),          intent(IN), value :: C_lbdd_CC_1
    type(C_PTR),          intent(IN), value :: C_CC_1
    type(C_PTR),          intent(IN), value :: C_lbdd_scratch_hydro_op1_auxC
    type(C_PTR),          intent(IN), value :: C_scratch_hydro_op1_auxC
    type(C_PTR),          intent(IN), value :: C_lbdd_scratch_hydro_op1_flX
    type(C_PTR),          intent(IN), value :: C_scratch_hydro_op1_flX
    type(C_PTR),          intent(IN), value :: C_lbdd_scratch_hydro_op1_flY
    type(C_PTR),          intent(IN), value :: C_scratch_hydro_op1_flY
    type(C_PTR),          intent(IN), value :: C_scratch_hydro_op1_flZ

    real             :: F_external_hydro_op1_dt
    integer          :: F_external_hydro_op1_eosMode
    real,    pointer :: F_tile_deltas(:)
    integer, pointer :: F_tile_hi(:)
    integer, pointer :: F_tile_interior(:, :)
    integer, pointer :: F_tile_lo(:)
    real,    pointer :: F_CC_1(:, :, :, :)
    real,    pointer :: F_scratch_hydro_op1_auxC(:, :, :)
    real,    pointer :: F_scratch_hydro_op1_flX(:, :, :, :)
    real,    pointer :: F_scratch_hydro_op1_flY(:, :, :, :)
    real,    pointer :: F_scratch_hydro_op1_flZ(:, :, :, :)
    integer, pointer :: F_lbdd_CC_1(:)
    integer, pointer :: F_lbdd_scratch_hydro_op1_auxC(:)
    integer, pointer :: F_lbdd_scratch_hydro_op1_flX(:)
    integer, pointer :: F_lbdd_scratch_hydro_op1_flY(:)

    F_external_hydro_op1_dt = REAL(MH_external_hydro_op1_dt)
    F_external_hydro_op1_eosMode = INT(MH_external_hydro_op1_eosMode)

    CALL C_F_POINTER(C_tile_deltas, F_tile_deltas, shape=[MILHOJA_MDIM])
    CALL C_F_POINTER(C_tile_hi, F_tile_hi, shape=[MILHOJA_MDIM])
    CALL C_F_POINTER(C_tile_interior, F_tile_interior, shape=[2, MILHOJA_MDIM])
    CALL C_F_POINTER(C_tile_lo, F_tile_lo, shape=[MILHOJA_MDIM])
    CALL C_F_POINTER(C_lbdd_CC_1, F_lbdd_CC_1, shape=[MILHOJA_MDIM+1])
    ! No need for variable masking in tile wrappers since the full data array is
    ! in the memory system and available for use.  Always include all variables?
    ! Makes sense, but how to get that number?
    CALL C_F_POINTER(C_CC_1, F_CC_1, shape=[8 + 2 * 1 * MILHOJA_K1D, 8 + 2 * 1 * MILHOJA_K2D, 1 + 2 * 1 * MILHOJA_K3D, 9 + 1 - 0])
    CALL C_F_POINTER(C_lbdd_scratch_hydro_op1_auxC, F_lbdd_scratch_hydro_op1_auxC, shape=[3])
    CALL C_F_POINTER(C_scratch_hydro_op1_auxC, F_scratch_hydro_op1_auxC, shape=[10, 10, 1])
    CALL C_F_POINTER(C_lbdd_scratch_hydro_op1_flX, F_lbdd_scratch_hydro_op1_flX, shape=[MILHOJA_MDIM+1])
    CALL C_F_POINTER(C_lbdd_scratch_hydro_op1_flY, F_lbdd_scratch_hydro_op1_flY, shape=[MILHOJA_MDIM+1])
    CALL C_F_POINTER(C_scratch_hydro_op1_flX, F_scratch_hydro_op1_flX, shape=[9, 8, 1, 5])
    CALL C_F_POINTER(C_scratch_hydro_op1_flY, F_scratch_hydro_op1_flY, shape=[8, 9, 1, 5])
    CALL C_F_POINTER(C_scratch_hydro_op1_flZ, F_scratch_hydro_op1_flZ, shape=[1, 1, 1, 1])

    !write(*,*) ""
    !write(*,*) "dt      = ", F_external_hydro_op1_dt
    !write(*,*) "eosMode = ", F_external_hydro_op1_eosMode
    !write(*,*) "deltas  = ", F_tile_deltas
    !write(*,*) "lo      = ", F_tile_interior(1, :)
    !write(*,*) "hi      = ", F_tile_interior(2, :)
    !write(*,*) "loU     = ", F_lbdd_CC_1
    !write(*,*) "loAuxC  = ", F_lbdd_scratch_hydro_op1_auxC
    !write(*,*) "loFl    = ", F_lbdd_scratch_hydro_op1_flX
    !write(*,*) ""

    CALL cpu_tf_hydro_Fortran( &
        F_external_hydro_op1_dt, &
        F_external_hydro_op1_eosMode, &
        F_tile_deltas, &
        F_tile_hi, &
        F_tile_interior, &
        F_tile_lo, &
        F_CC_1, &
        F_scratch_hydro_op1_auxC, &
        F_scratch_hydro_op1_flX, &
        F_scratch_hydro_op1_flY, &
        F_scratch_hydro_op1_flZ, &
        F_lbdd_CC_1, &
        F_lbdd_scratch_hydro_op1_auxC, &
        F_lbdd_scratch_hydro_op1_flX, &
        F_lbdd_scratch_hydro_op1_flY)
end subroutine cpu_tf_hydro_C2F
