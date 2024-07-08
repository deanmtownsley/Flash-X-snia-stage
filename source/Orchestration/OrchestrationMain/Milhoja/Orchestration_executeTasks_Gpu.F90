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
#include "Milhoja.h"

#ifdef ORCHESTRATION_USE_DATAPACKETS
!> @ingroup OrchestrationMilhoja
!! @stubref{Orchestration_executeTasks_Gpu}
!!
!! @brief Concrete implementation of Orchestration_executeTasks_Gpu
subroutine Orchestration_executeTasks_Gpu(MH_taskFunction,     &
                                          nDistributorThreads, &
                                          nThreads,            &
                                          nTilesPerPacket,     &
                                          MH_packet_Cptr)
    use iso_c_binding,       ONLY : C_PTR

    use milhoja_types_mod,   ONLY : MILHOJA_INT
    use milhoja_runtime_mod, ONLY : milhoja_runtime_taskFunction
#ifdef RUNTIME_CAN_USE_TILEITER
    use milhoja_runtime_mod, ONLY : milhoja_runtime_executeTasks_Gpu
#endif

    use Orchestration_interface, ONLY : Orchestration_checkInternalError

    implicit none

    procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
    integer,                                intent(IN) :: nDistributorThreads
    integer,                                intent(IN) :: nThreads
    integer,                                intent(IN) :: nTilesPerPacket
    type(C_PTR),                            intent(IN) :: MH_packet_CPtr

    integer(MILHOJA_INT) :: MH_nDistributorThreads
    integer(MILHOJA_INT) :: MH_nThreads
    integer(MILHOJA_INT) :: MH_nTilesPerPacket
    integer(MILHOJA_INT) :: MH_ierr

    MH_nDistributorThreads = INT(nDistributorThreads, kind=MILHOJA_INT)
    MH_nThreads            = INT(nThreads,            kind=MILHOJA_INT)
    MH_nTilesPerPacket     = INT(nTilesPerPacket,     kind=MILHOJA_INT)

#ifdef RUNTIME_CAN_USE_TILEITER
    CALL milhoja_runtime_executeTasks_Gpu(MH_taskFunction,        &
                                          MH_nDistributorThreads, &
                                          MH_nThreads,            &
                                          MH_nTilesPerPacket,     &
                                          MH_packet_Cptr,         &
                                          MH_ierr)
    CALL Orchestration_checkInternalError("Orchestration_executeTasks_Gpu", MH_ierr)
#else
    CALL Driver_abort("Orchestration_executeTasks_Gpu: milhoja_runtime_executeTasks_Gpu disabled")
#endif
end subroutine Orchestration_executeTasks_Gpu
#endif

