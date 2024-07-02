!> @copyright Copyright 2024 UChicago Argonne, LLC and contributors
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
!! interface.
module Orchestration_interface

    use Orchestration_interfaceTypeDecl, ONLY: Orchestration_tileCInfo_t, &
                                               MILHOJA_INT
    use iso_c_binding, ONLY : C_PTR
    implicit none

    interface
        subroutine Orchestration_init()
            implicit none
        end subroutine Orchestration_init

        subroutine Orchestration_finalize()
            implicit none
        end subroutine Orchestration_finalize

        subroutine Orchestration_checkInternalError(routineName, MH_errorCode)
            import
            implicit none
            character(LEN=*),     intent(IN) :: routineName
            integer(MILHOJA_INT), intent(IN) :: MH_errorCode
        end subroutine Orchestration_checkInternalError
    end interface

    interface
        subroutine Orchestration_setupPipelineForCpuTasks(MH_taskFunction, nThreads)
            use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
            integer,                                intent(IN) :: nThreads
        end subroutine Orchestration_setupPipelineForCpuTasks

        subroutine Orchestration_setupPipelineForGpuTasks(MH_taskFunction, &
                                          nThreads,            &
                                          nTilesPerPacket,     &
                                          MH_packet_Cptr)
            use iso_c_binding, ONLY : C_PTR
            use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
            integer,                                intent(IN) :: nThreads
            integer,                                intent(IN) :: nTilesPerPacket
            type(C_PTR),                            intent(IN) :: MH_packet_CPtr
        end subroutine Orchestration_setupPipelineForGpuTasks

        ! Note that the following is for both the CPU/GPU and Split
        ! CPU/CPU bundles, depending on presence of the
        ! nTilesPerCpuTurn argument. - KW
        subroutine Orchestration_setupPipelineForCpuGpuTasks(MH_pktTaskFunction, &
             MH_tileTaskFunction, &
             nThreads,            &
             nTilesPerPacket,     &
             MH_pktProto_Cptr,    &
             nTilesPerCpuTurn)
            use iso_c_binding, ONLY : C_PTR
            use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_pktTaskFunction
            procedure(milhoja_runtime_taskFunction)            :: MH_tileTaskFunction
            integer,                                intent(IN) :: nThreads
            integer,                                intent(IN) :: nTilesPerPacket
            type(C_PTR),                            intent(IN) :: MH_pktProto_CPtr
            integer,OPTIONAL,                       intent(IN) :: nTilesPerCpuTurn
        end subroutine Orchestration_setupPipelineForCpuGpuTasks

        subroutine Orchestration_setupPipelineForExtGpuTasks(MH_taskFunction, &
                                          MH_postTaskFunction,                &
                                          nThreads,            &
                                          nTilesPerPacket,     &
                                          MH_packet_Cptr,     &
                                          MH_postProto_Cptr)
            use iso_c_binding, ONLY : C_PTR
            use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
            procedure(milhoja_runtime_taskFunction)            :: MH_postTaskFunction
            integer,                                intent(IN) :: nThreads
            integer,                                intent(IN) :: nTilesPerPacket
            type(C_PTR),                            intent(IN) :: MH_packet_CPtr
            type(C_PTR),                            intent(IN) :: MH_postProto_CPtr
        end subroutine Orchestration_setupPipelineForExtGpuTasks
    end interface

    interface
        subroutine Orchestration_executeTasks_Cpu(MH_taskFunction, &
                                                  prototype_Cptr, nThreads)
            use iso_c_binding, ONLY : C_PTR
            use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
            type(C_PTR),                            intent(IN) :: prototype_Cptr
            integer,                                intent(IN) :: nThreads
        end subroutine Orchestration_executeTasks_Cpu


        subroutine Orchestration_executeTasks_Gpu(MH_taskFunction,     &
                                                  nDistributorThreads, &
                                                  nThreads,            &
                                                  nTilesPerPacket,     &
                                                  MH_packet_Cptr)
            use iso_c_binding,       ONLY : C_PTR
            use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
            integer,                                intent(IN) :: nDistributorThreads
            integer,                                intent(IN) :: nThreads
            integer,                                intent(IN) :: nTilesPerPacket
            type(C_PTR),                            intent(IN) :: MH_packet_Cptr
        end subroutine Orchestration_executeTasks_Gpu
        subroutine Orchestration_executeTasks_CpuGpu(MH_pktTaskFunction, &
                                                  MH_tileTaskFunction, &
                                                  nDistributorThreads, &
                                                  nThreads,            &
                                                  nTilesPerPacket,     &
                                                  MH_packet_Cptr,     &
                                                  MH_tileProto_Cptr)
            use iso_c_binding,       ONLY : C_PTR
            use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_pktTaskFunction
            procedure(milhoja_runtime_taskFunction)            :: MH_tileTaskFunction
            integer,                                intent(IN) :: nDistributorThreads
            integer,                                intent(IN) :: nThreads
            integer,                                intent(IN) :: nTilesPerPacket
            type(C_PTR),                            intent(IN) :: MH_packet_Cptr
            type(C_PTR),                            intent(IN) :: MH_tileProto_Cptr
        end subroutine Orchestration_executeTasks_CpuGpu
        subroutine Orchestration_executeTasks_CpuGpuSplit(MH_pktTaskFunction, &
                                                  MH_tileTaskFunction, &
                                                  nDistributorThreads, &
                                                  nThreads,            &
                                                  nTilesPerPacket,     &
                                                  nTilesPerCpuTurn,    &
                                                  MH_packet_Cptr,     &
                                                  MH_tileProto_Cptr)
            use iso_c_binding,       ONLY : C_PTR
            use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_pktTaskFunction
            procedure(milhoja_runtime_taskFunction)            :: MH_tileTaskFunction
            integer,                                intent(IN) :: nDistributorThreads
            integer,                                intent(IN) :: nThreads
            integer,                                intent(IN) :: nTilesPerPacket
            integer,                                intent(IN) :: nTilesPerCpuTurn
            type(C_PTR),                            intent(IN) :: MH_packet_Cptr
            type(C_PTR),                            intent(IN) :: MH_tileProto_Cptr
        end subroutine Orchestration_executeTasks_CpuGpuSplit
        subroutine Orchestration_executeTasks_extGpu(MH_taskFunction,  &
                                                  MH_postTaskFunction, &
                                                  nDistributorThreads, &
                                                  nThreads,            &
                                                  nTilesPerPacket,     &
                                                  MH_packet_Cptr,     &
                                                  MH_postProto_Cptr)
            use iso_c_binding,       ONLY : C_PTR
            use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction
            implicit none
            procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
            procedure(milhoja_runtime_taskFunction)            :: MH_postTaskFunction
            integer,                                intent(IN) :: nDistributorThreads
            integer,                                intent(IN) :: nThreads
            integer,                                intent(IN) :: nTilesPerPacket
            type(C_PTR),                            intent(IN) :: MH_packet_Cptr
            type(C_PTR),                            intent(IN) :: MH_postProto_Cptr
        end subroutine Orchestration_executeTasks_extGpu
    end interface


    !Separate specific interfaces - push to pipeline
    interface
        !The first one here is for a packet-less CPU-only pipeline
        subroutine Orchestration_pushTileToCpuPipeline(prototype_Cptr, nThreads, &
                                                    tileCInfo)
            import
            implicit none
            type(C_PTR),                            intent(IN) :: prototype_Cptr
            integer,                                intent(IN) :: nThreads
            type(Orchestration_tileCInfo_t),target, intent(IN) :: tileCInfo
        end subroutine Orchestration_pushTileToCpuPipeline
        subroutine Orchestration_pushTileToGpuPipeline(prototype_Cptr, nThreads, &
                                                    tileCInfo)
            import
            implicit none
            type(C_PTR),                            intent(IN) :: prototype_Cptr
            integer,                                intent(IN) :: nThreads
            type(Orchestration_tileCInfo_t),target, intent(IN) :: tileCInfo
        end subroutine Orchestration_pushTileToGpuPipeline
        subroutine Orchestration_pushTileToCpuGpuPipeline(pktProto_Cptr, tileProto_Cptr, &
                                                    nThreads, &
                                                    tileCInfo)
            import
            implicit none
            type(C_PTR),                            intent(IN) :: pktProto_Cptr
            type(C_PTR),                            intent(IN) :: tileProto_Cptr
            integer,                                intent(IN) :: nThreads
            type(Orchestration_tileCInfo_t),target, intent(IN) :: tileCInfo
        end subroutine Orchestration_pushTileToCpuGpuPipeline
        subroutine Orchestration_pushTileToCpuGpuSplitPipeline(pktProto_Cptr, tileProto_Cptr, &
                                                    nThreads, &
                                                    tileCInfo)
            import
            implicit none
            type(C_PTR),                            intent(IN) :: pktProto_Cptr
            type(C_PTR),                            intent(IN) :: tileProto_Cptr
            integer,                                intent(IN) :: nThreads
            type(Orchestration_tileCInfo_t),target, intent(IN) :: tileCInfo
        end subroutine Orchestration_pushTileToCpuGpuSplitPipeline
        subroutine Orchestration_pushTileToExtGpuPipeline(prototype_Cptr, nThreads, &
                                                    tileCInfo)
            import
            implicit none
            type(C_PTR),                            intent(IN) :: prototype_Cptr
            integer,                                intent(IN) :: nThreads
            type(Orchestration_tileCInfo_t),target, intent(IN) :: tileCInfo
        end subroutine Orchestration_pushTileToExtGpuPipeline
    end interface

    ! Generic interface - tear down pipeline
    interface Orchestration_teardownPipeline
        subroutine Orchestration_teardownPipelineForCpuTasks(nThreads)
            implicit none
            integer,                              intent(IN) :: nThreads
        end subroutine Orchestration_teardownPipelineForCpuTasks
        subroutine Orchestration_teardownPipelineForGpuTasks(nThreads, nTilesPerPacket)
            implicit none
            integer,                              intent(IN) :: nThreads
            integer,                              intent(IN) :: nTilesPerPacket
        end subroutine Orchestration_teardownPipelineForGpuTasks
    end interface Orchestration_teardownPipeline
    interface
        subroutine Orchestration_teardownPipelineForCpuGpuTasks(nThreads, nTilesPerPacket)
            implicit none
            integer,                              intent(IN) :: nThreads
            integer,                              intent(IN) :: nTilesPerPacket
        end subroutine Orchestration_teardownPipelineForCpuGpuTasks
    end interface
    interface
        subroutine Orchestration_teardownPipelineForCpuGpuSplitTasks(nThreads, nTilesPerPacket)
            implicit none
            integer,                              intent(IN) :: nThreads
            integer,                              intent(IN) :: nTilesPerPacket
        end subroutine Orchestration_teardownPipelineForCpuGpuSplitTasks
    end interface
    interface
        subroutine Orchestration_teardownPipelineForExtGpuTasks(nThreads, nTilesPerPacket)
            implicit none
            integer,                              intent(IN) :: nThreads
            integer,                              intent(IN) :: nTilesPerPacket
        end subroutine Orchestration_teardownPipelineForExtGpuTasks
    end interface

end module Orchestration_interface
! Local Variables:
! f90-program-indent: 4
! f90-do-indent: 4
! f90-type-indent: 4
! indent-tabs-mode: nil
! End:
