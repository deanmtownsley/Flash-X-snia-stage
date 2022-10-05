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

#include "Simulation.h"

!> @ingroup Orchestration
!!
!! @brief Public interface of Orchestration unit
!!
!! @details
!! A standard Flash-X Fortran module that encapsulates the interface declarations
!! of all routine's in the Orchestration unit that are part of this unit's public
!! interface.  Refer to the doxygen documentation for this unit for a discussion of
!! why the interface changes based on setup information.
module Orchestration_interface

    implicit none

    interface
        subroutine Orchestration_init()
            implicit none
        end subroutine Orchestration_init

        subroutine Orchestration_finalize()
            implicit none
        end subroutine Orchestration_finalize

#ifdef USE_MILHOJA_RUNTIME
        subroutine Orchestration_checkInternalError(routineName, MH_errorCode)
            use milhoja_types_mod, ONLY : MILHOJA_INT
            implicit none
            character(LEN=*),     intent(IN) :: routineName
            integer(MILHOJA_INT), intent(IN) :: MH_errorCode
        end subroutine Orchestration_checkInternalError

        subroutine Orchestration_executeTasks_Cpu(MH_taskFunction, nThreads)
            use milhoja_runtime_mod, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
            integer,                                intent(IN) :: nThreads
        end subroutine Orchestration_executeTasks_Cpu

#ifdef ORCHESTRATION_USE_GPUS
        subroutine Orchestration_executeTasks_Gpu(MH_taskFunction,     &
                                                  nDistributorThreads, &
                                                  nThreads,            &
                                                  nTilesPerPacket,     &
                                                  MH_packet_Cptr)
            use iso_c_binding,       ONLY : C_PTR
            use milhoja_runtime_mod, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
            integer,                                intent(IN) :: nDistributorThreads
            integer,                                intent(IN) :: nThreads
            integer,                                intent(IN) :: nTilesPerPacket
            type(C_PTR),                            intent(IN) :: MH_packet_Cptr
        end subroutine Orchestration_executeTasks_Gpu
#endif
#endif
    end interface

end module Orchestration_interface

