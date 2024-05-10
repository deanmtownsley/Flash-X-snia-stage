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
!! @brief Concrete implementation of Orchestration_setupPipelineForCpuGpuTasks
#include "Milhoja.h"
subroutine Orchestration_setupPipelineForCpuGpuTasks(MH_pktTaskFunction, &
                                          MH_tileTaskFunction,     &
                                          nThreads,            &
                                          nTilesPerPacket,     &
                                          MH_pktProto_Cptr)
    use iso_c_binding, ONLY : C_PTR

    use milhoja_types_mod,   ONLY : MILHOJA_INT
    use milhoja_runtime_mod, ONLY : milhoja_runtime_taskFunction
#if ! defined(RUNTIME_MUST_USE_TILEITER) && defined(RUNTIME_SUPPORT_DATAPACKETS)
    use milhoja_runtime_mod, ONLY : milhoja_runtime_setupPipelineForCpuGpuTasks
#endif

    use Driver_interface,        ONLY : Driver_abort
    use Orchestration_interface, ONLY : Orchestration_checkInternalError

    implicit none

    procedure(milhoja_runtime_taskFunction)            :: MH_pktTaskFunction
    procedure(milhoja_runtime_taskFunction)            :: MH_tileTaskFunction
    integer,                                intent(IN) :: nThreads
    integer,                                intent(IN) :: nTilesPerPacket
    type(C_PTR),                            intent(IN) :: MH_pktProto_CPtr

    integer(MILHOJA_INT) :: MH_nThreads
    integer(MILHOJA_INT) :: MH_nTilesPerPacket
    integer(MILHOJA_INT) :: MH_ierr

    MH_nThreads = INT(nThreads, kind=MILHOJA_INT)
    MH_nTilesPerPacket     = INT(nTilesPerPacket,     kind=MILHOJA_INT)

#if ! defined(RUNTIME_MUST_USE_TILEITER) && defined(RUNTIME_SUPPORT_DATAPACKETS)
    CALL milhoja_runtime_setupPipelineForCpuGpuTasks(MH_pktTaskFunction, &
                                          MH_tileTaskFunction,        &
                                          MH_nThreads,            &
                                          MH_nTilesPerPacket,     &
                                          MH_pktProto_Cptr,         &
                                          MH_ierr)
    CALL Orchestration_checkInternalError("Orchestration_setupPipelineForCpuGpuTasks", MH_ierr)
#else
    CALL Driver_abort("Orchestration_setupPipelineForCpuGpuTasks: milhoja_runtime_setupPipelineForCpuGpuTasks disabled")
#endif
end subroutine Orchestration_setupPipelineForCpuGpuTasks

