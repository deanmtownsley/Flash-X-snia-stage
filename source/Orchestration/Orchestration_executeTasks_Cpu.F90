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

#ifdef FLASHX_ORCHESTRATION_MILHOJA
!> @ingroup Orchestration
!! @anchor Orchestration_executeTasks_Cpu_stub
!!
!! @brief Execute a task function with the CPU-only thread team configuration
!!
!! @details
!! Use the runtime to execute the given task function using the CPU-only thread
!! team configuration.  Upon termination, the task function will have been applied
!! by the CPU to all leaf blocks using the given number of threads.  The order
!! in which the task function is applied to blocks is determined at runtime and
!! should therefore be considered as arbitrary.  It is a logical error to call this
!! routine if the runtime is not enabled.
!!
!! @todo Presently, the names of executeTasks routines map onto the associated
!! thread team configuration.  We need a practical and scalable scheme to
!! express and manage this mapping.
!! @todo This interface is presently restricted to leaf blocks only because the
!! runtime implementation is limited in this way.  Once the runtime functions with
!! AMR, this interface will need to include a specification of which blocks to
!! apply the TF to as well as if tiling should be used.
!! @todo For logging purposes, this interface should accept a name for the task
!! function.
!!
!! @param MH_taskFunction  The task function to be executed by the single
!!                         thread team.  It is assumed that this function
!!                         was written to run on CPUs.
!! @param nThreads         The number of threads to activate in team
subroutine Orchestration_executeTasks_Cpu(MH_taskFunction, nThreads)
    use milhoja_runtime_mod, ONLY : milhoja_runtime_taskFunction
    use Driver_interface,    ONLY : Driver_abort

    implicit none

    procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
    integer,                                intent(IN) :: nThreads

    CALL Driver_abort("[Orchestration_executeTasks_Cpu] Runtime not enabled")
end subroutine Orchestration_executeTasks_Cpu
#endif

