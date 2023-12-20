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

!> @details
!! A module that provides Fortran interfaces to all C++ functions written for
!! working with this particular data item.  This module assumes that all these
!! functions were written with C-compatible linkage.
module Tile_cpu_tf_hydro_C2F_mod
    implicit none
    private

    public :: new_hydro_advance_wrapper_C
    public :: delete_hydro_advance_wrapper_C
    public :: acquire_scratch_wrapper_C
    public :: release_scratch_wrapper_C

    !!!!!----- INTERFACES TO C-LINKAGE C++ FUNCTIONS FOR TIME ADVANCE UNIT
    interface
        !> To be used by TimeAdvance to create a concrete tile wrapper that the
        !! Orchestation unit can use to blindly clone the same type of wrapper.
        function new_hydro_advance_wrapper_C(C_external_hydro_op1_dt, &
                                             C_external_hydro_op1_eosMode, &
                                             C_wrapper) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT, &
                                          MILHOJA_REAL
            real(MILHOJA_REAL),   intent(IN), value :: C_external_hydro_op1_dt
            integer(MILHOJA_INT), intent(IN), value :: C_external_hydro_op1_eosMode
            type(C_PTR),          intent(IN)        :: C_wrapper
            integer(MILHOJA_INT)                    :: C_ierr
        end function new_hydro_advance_wrapper_C

        !> To be used by TimeAdvance to free tile wrapper resources.
        function delete_hydro_advance_wrapper_C(C_wrapper) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            type(C_PTR),         intent(IN), value :: C_wrapper
            integer(MILHOJA_INT)                   :: C_ierr
        end function delete_hydro_advance_wrapper_C

        !> To be used by TimeAdvance to acquire thread-private scratch before
        !! the prototype TileWrapper object is used
        function acquire_scratch_wrapper_C() result(C_ierr) bind(c)
            use milhoja_types_mod, ONLY : MILHOJA_INT
            integer(MILHOJA_INT) :: C_ierr
        end function acquire_scratch_wrapper_C

        !> To be used by TimeAdvance to release thread-private scratch after the
        !! prototype TileWrapper object has been used
        function release_scratch_wrapper_C() result(C_ierr) bind(c)
            use milhoja_types_mod, ONLY : MILHOJA_INT
            integer(MILHOJA_INT) :: C_ierr
        end function release_scratch_wrapper_C
    end interface

end module Tile_cpu_tf_hydro_C2F_mod
