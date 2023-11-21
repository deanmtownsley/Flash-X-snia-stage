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
module DataPacket_gpu_tf_hydro_C2F_mod
    implicit none
    private

    public :: instantiate_gpu_tf_hydro_packet_C
    public :: delete_gpu_tf_hydro_packet_C
    public :: release_gpu_tf_hydro_extra_queue_C

    !!!!!----- INTERFACES TO C-LINKAGE C++ FUNCTIONS FOR TIME ADVANCE UNIT
    interface
        !> To be used by TimeAdvance to create a concrete data packet that the
        !! Orchestration unit can use to blindly clone the same type of packet.
        function instantiate_gpu_tf_hydro_packet_C(C_dt, C_packet) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT, &
                                          MILHOJA_REAL
            real(MILHOJA_REAL),  intent(IN), value :: C_dt
            type(C_PTR),         intent(IN)        :: C_packet
            integer(MILHOJA_INT)                   :: C_ierr
        end function instantiate_gpu_tf_hydro_packet_C

        !> To be used by TimeAdvance to free packet resources.
        function delete_gpu_tf_hydro_packet_C(C_packet) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            type(C_PTR),         intent(IN), value :: C_packet
            integer(MILHOJA_INT)                   :: C_ierr
        end function delete_gpu_tf_hydro_packet_C
    end interface

    !!!!!----- INTERFACES TO C-LINKAGE C++ FUNCTIONS FOR TASK FUNCTION
    interface
        !> To be used by task function to release extra OpenACC asynchronous
        !! queue resources
        function release_gpu_tf_hydro_extra_queue_C(C_packet, C_id) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            type(C_PTR),          intent(IN), value :: C_packet
            integer(MILHOJA_INT), intent(IN), value :: C_id
            integer(MILHOJA_INT)                    :: C_ierr
        end function release_gpu_tf_hydro_extra_queue_C
    end interface

end module DataPacket_gpu_tf_hydro_C2F_mod
