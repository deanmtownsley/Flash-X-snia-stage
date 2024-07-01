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
!! @anchor Orchestration_executeTasks_Cpu_stub
!!
!! @brief Execute a task function with the CPU-only thread team configuration
!!
!! @details
!! Use the runtime to execute the given task function using the CPU-only thread
!! team configuration.  Upon termination, the task function will have been applied
!! by the CPU to all leaf blocks using the given number of threads.  The order
!! in which the task function is applied to blocks is determined at runtime and
!! should therefore be considered as arbitrary.
!!
!! @note
!! An actual, non-stub implementation of this interface is only available if (1)
!! Flash-X is linked with a Milhoja library that was configured to support
!! wrapped tiles and (2) the linked Milhoja library was also configured to support
!! execute-style orchestration calls.
!! The only existing implementation additionally requires that Flash-X uses the
!! "Milhoja" Grid implementation and that the Milhoja library fully provides the
!! necessary grid backend; currently known variants of the Milhoja library do
!! this by using the AMReX library as their grid backend.
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
!! @param prototype_Cptr   Pointer to a prototype tile wrapper object to be
!!                         used to create new wrapped tiles.
!! @param nThreads         The number of threads to activate in team
subroutine Orchestration_executeTasks_Cpu(MH_taskFunction, &
                                          prototype_Cptr, nThreads)
    use iso_c_binding, ONLY : C_PTR

    use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction
    use Driver_interface,    ONLY : Driver_abort

    implicit none

    procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
    type(C_PTR),                            intent(IN) :: prototype_Cptr
    integer,                                intent(IN) :: nThreads

#ifdef FLASHX_ORCHESTRATION_MILHOJA
    CALL Driver_abort("[Orchestration_executeTasks_Cpu] Runtime not enabled")
#endif
end subroutine Orchestration_executeTasks_Cpu
