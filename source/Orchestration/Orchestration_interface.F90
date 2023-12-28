!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
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

    use Orchestration_interfaceTypeDecl, ONLY: Orchestration_tileCInfo_t
    implicit none

    interface
        subroutine Orchestration_init()
            implicit none
        end subroutine Orchestration_init

        subroutine Orchestration_finalize()
            implicit none
        end subroutine Orchestration_finalize

#ifdef FLASHX_ORCHESTRATION_MILHOJA
        subroutine Orchestration_checkInternalError(routineName, MH_errorCode)
            use milhoja_types_mod, ONLY : MILHOJA_INT
            implicit none
            character(LEN=*),     intent(IN) :: routineName
            integer(MILHOJA_INT), intent(IN) :: MH_errorCode
        end subroutine Orchestration_checkInternalError

#ifndef RUNTIME_USES_TILEITER
        subroutine Orchestration_setupPipelineForCpuTasks(MH_taskFunction, &
                                                  prototype_Cptr, nThreads)
            use iso_c_binding, ONLY : C_PTR
            use milhoja_runtime_mod, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
            type(C_PTR),                            intent(IN) :: prototype_Cptr
            integer,                                intent(IN) :: nThreads
        end subroutine Orchestration_setupPipelineForCpuTasks

        subroutine Orchestration_teardownPipelineForCpuTasks(MH_taskFunction, &
                                                  prototype_Cptr, nThreads)
            use iso_c_binding, ONLY : C_PTR
            use milhoja_runtime_mod, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
            type(C_PTR),                            intent(IN) :: prototype_Cptr
            integer,                                intent(IN) :: nThreads
        end subroutine Orchestration_teardownPipelineForCpuTasks

        subroutine Orchestration_pushTileToPipeline(MH_taskFunction, &
                                                    prototype_Cptr, nThreads, &
                                                    tileCInfo)
            use iso_c_binding, ONLY : C_PTR
            use milhoja_runtime_mod, ONLY : milhoja_runtime_taskFunction
            import
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
            type(C_PTR),                            intent(IN) :: prototype_Cptr
            integer,                                intent(IN) :: nThreads
            type(Orchestration_tileCInfo_t),        intent(IN) :: tileCInfo
        end subroutine Orchestration_pushTileToPipeline

#else
        subroutine Orchestration_executeTasks_Cpu(MH_taskFunction, &
                                                  prototype_Cptr, nThreads)
            use iso_c_binding, ONLY : C_PTR
            use milhoja_runtime_mod, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
            type(C_PTR),                            intent(IN) :: prototype_Cptr
            integer,                                intent(IN) :: nThreads
        end subroutine Orchestration_executeTasks_Cpu
#endif

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
! Local Variables:
! f90-program-indent: 4
! f90-do-indent: 4
! f90-type-indent: 4
! f90-associate-indent: 9
! f90-if-indent: 40
! indent-tabs-mode: nil
! End:
