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
!! @stubref{Orchestration_teardownPipelineForCpuGpuSplitTasks}
!!
!! @brief Concrete implementation of Orchestration_teardownPipelineForCpuGpuSplitTasks
#include "Milhoja.h"
subroutine Orchestration_teardownPipelineForCpuGpuSplitTasks(nThreads, nTilesPerPacket)
    use milhoja_types_mod,   ONLY : MILHOJA_INT
#if ! defined(RUNTIME_MUST_USE_TILEITER) && defined(RUNTIME_SUPPORT_DATAPACKETS)
    use milhoja_runtime_mod, ONLY : milhoja_runtime_teardownPipelineForCpuGpuSplitTasks
#endif

    use Driver_interface,        ONLY : Driver_abort
    use Orchestration_interface, ONLY : Orchestration_checkInternalError

    implicit none

    integer,                                intent(IN) :: nThreads
    integer,                                intent(IN) :: nTilesPerPacket

    integer(MILHOJA_INT) :: MH_nThreads
    integer(MILHOJA_INT) :: MH_nTilesPerPacket
    integer(MILHOJA_INT) :: MH_ierr

    MH_nThreads = INT(nThreads, kind=MILHOJA_INT)
    MH_nTilesPerPacket = INT(nTilesPerPacket, kind=MILHOJA_INT)

#if ! defined(RUNTIME_MUST_USE_TILEITER) && defined(RUNTIME_SUPPORT_DATAPACKETS)
    CALL milhoja_runtime_teardownPipelineForCpuGpuSplitTasks(MH_nThreads, MH_nTilesPerPacket, MH_ierr)
    CALL Orchestration_checkInternalError("Orchestration_teardownPipelineForCpuGpuSplitTasks", MH_ierr)
#else
    CALL Driver_abort("Orchestration_teardownPipelineForCpuGpuSplitTasks:&
         & milhoja_runtime_teardownPipelineForCpuGpuSplitTasks disabled")
#endif
end subroutine Orchestration_teardownPipelineForCpuGpuSplitTasks

