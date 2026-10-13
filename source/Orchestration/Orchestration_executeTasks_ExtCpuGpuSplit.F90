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

!> @ingroup OrchestrationMilhoja
!! @stubref{Orchestration_executeTasks_ExtCpuGpuSplit}
!!
!! @brief Concrete implementation of Orchestration_executeTasks_ExtCpuGpuSplit
subroutine Orchestration_executeTasks_ExtCpuGpuSplit(MH_tileTaskFunction,     &
                                                     MH_pktTaskFunction,      &
                                                     MH_postTileTaskFunction, &
                                                     nDistributorThreads,     &
                                                     nThreads,                &
                                                     nTilesPerPacket,         &
                                                     nTilesPerCpuTurn,        &
                                                     MH_packet_Cptr,          &
                                                     MH_tileProto_Cptr,       &
                                                     MH_postTileProto_Cptr)
    use iso_c_binding, ONLY : C_PTR
    use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction

    implicit none

    procedure(milhoja_runtime_taskFunction)            :: MH_tileTaskFunction
    procedure(milhoja_runtime_taskFunction)            :: MH_pktTaskFunction
    procedure(milhoja_runtime_taskFunction)            :: MH_postTileTaskFunction
    integer,                                intent(IN) :: nDistributorThreads
    integer,                                intent(IN) :: nThreads
    integer,                                intent(IN) :: nTilesPerPacket
    integer,                                intent(IN) :: nTilesPerCpuTurn
    type(C_PTR),                            intent(IN) :: MH_packet_Cptr
    type(C_PTR),                            intent(IN) :: MH_tileProto_Cptr
    type(C_PTR),                            intent(IN) :: MH_postTileProto_Cptr

end subroutine Orchestration_executeTasks_ExtCpuGpuSplit

