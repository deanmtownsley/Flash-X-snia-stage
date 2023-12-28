!! @copyright Copyright 2023 UChicago Argonne, LLC and contributors
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
!! @stubref{Orchestration_setupPipelineForCpuTasks}
!!
!! @brief Concrete implementation of Orchestration_setupPipelineForCpuTasks
subroutine Orchestration_setupPipelineForCpuTasks(MH_taskFunction, &
                                          prototype_Cptr, nThreads)
    use iso_c_binding, ONLY : C_PTR

    use milhoja_types_mod,   ONLY : MILHOJA_INT
    use milhoja_runtime_mod, ONLY : milhoja_runtime_taskFunction, &
                                    milhoja_runtime_setupPipelineForCpuTasks

    use Orchestration_interface, ONLY : Orchestration_checkInternalError

    implicit none

    procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
    type(C_PTR),                            intent(IN) :: prototype_Cptr
    integer,                                intent(IN) :: nThreads

    integer(MILHOJA_INT) :: MH_nThreads
    integer(MILHOJA_INT) :: MH_ierr

    MH_nThreads = INT(nThreads, kind=MILHOJA_INT)

    CALL milhoja_runtime_setupPipelineForCpuTasks(MH_taskFunction, prototype_Cptr, &
                                          MH_nThreads, MH_ierr)
    CALL Orchestration_checkInternalError("Orchestration_setupPipelineForCpuTasks", MH_ierr)
end subroutine Orchestration_setupPipelineForCpuTasks

