!! @copyright Copyright 2024 UChicago Argonne, LLC and contributors
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

!> @ingroup OrchestrationMilhoja
!! @stubref{Orchestration_setupPipelineForCpuGpuTasks}
!!
!! @brief Stub implementation of Orchestration_setupPipelineForCpuGpuTasks
#include "Milhoja.h"
subroutine Orchestration_setupPipelineForCpuGpuTasks(MH_pktTaskFunction, &
                                          MH_tileTaskFunction,     &
                                          nThreads,            &
                                          nTilesPerPacket,     &
                                          MH_pktProto_Cptr)
    use iso_c_binding, ONLY : C_PTR

    use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction

    implicit none

    procedure(milhoja_runtime_taskFunction)            :: MH_pktTaskFunction
    procedure(milhoja_runtime_taskFunction)            :: MH_tileTaskFunction
    integer,                                intent(IN) :: nThreads
    integer,                                intent(IN) :: nTilesPerPacket
    type(C_PTR),                            intent(IN) :: MH_pktProto_CPtr
end subroutine Orchestration_setupPipelineForCpuGpuTasks

