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

#if defined(FLASHX_ORCHESTRATION_MILHOJA) && defined(ORCHESTRATION_USE_GPUS)
!> @ingroup Orchestration
!! @anchor Orchestration_executeTasks_Gpu_stub
!!
!! @brief Execute a task function with the GPU-only thread team configuration
!!
!! @details
!! Use the runtime to execute the given task function with the GPU-only thread team
!! configuration.  Upon termination, the task function will have been applied to
!! all leaf blocks by the GPUs assigned to each MPI process' runtime instance and in
!! accord with the given runtime configuration values.  The order in which the
!! task function is applied to blocks is determined at runtime and should therefore
!! be considered as arbitrary.
!!
!! @note
!! This subroutine is only available in the Orchestration unit's public interface
!! if the Milhoja Orchestration implementation is in use and a GPU-compatible
!! runtime backend has been specified.  Therefore, it effectively has no stub
!! implementation.  Refer to the Orchestration unit's documentation for more
!! information.
!!
!! @todo This interface is presently restricted to leaf blocks only because the
!! runtime implementation is limited in this way.  Once the runtime functions with
!! AMR, this interface will need to include a specification of which blocks to
!! apply the TF to as well as if tiling should be used.
!! @todo For logging purposes, this interface should accept a name for the task
!! function.
!!
!! @param MH_taskFunction       The task function to be executed by the single
!!                              thread team.  It is assumed that this function
!!                              was written to launch kernels on GPUs
!! @param nDistributorThreads   The number of distributor threads to use
!! @param nThreads              The number of threads to activate in team
!! @param nTilesPerPacket       The maximum number of tiles allowed in each
!!                              packet
!! @param MH_packet_Cptr        Pointer to a prototype data packet to be
!!                              used to create new packets
subroutine Orchestration_executeTasks_Gpu(MH_taskFunction,     &
                                          nDistributorThreads, &
                                          nThreads,            &
                                          nTilesPerPacket,     &
                                          MH_packet_Cptr)
    use iso_c_binding,       ONLY : C_PTR

    use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction
    use Driver_interface,    ONLY : Driver_abort

    implicit none

    procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
    integer,                                intent(IN) :: nDistributorThreads
    integer,                                intent(IN) :: nThreads
    integer,                                intent(IN) :: nTilesPerPacket
    type(C_PTR),                            intent(IN) :: MH_packet_Cptr

    CALL Driver_abort("[Orchestration_executeTasks_Gpu] Runtime not enabled or GPUs not available")
end subroutine Orchestration_executeTasks_Gpu
#endif

