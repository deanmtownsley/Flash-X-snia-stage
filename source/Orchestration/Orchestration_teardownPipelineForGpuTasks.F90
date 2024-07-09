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
!! @stubref{Orchestration_teardownPipelineForGpuTasks}
!!
!! @brief Stub implementation of Orchestration_teardownPipelineForGpuTasks
!!
!! @details
!! Tear down the thread team bundle of the orchestration runtime that was
!! set up with Orchestration_setupPipelineForGpuTasks.
!! If necessary, queues of the bundle will be flushed and drained, and
!! pending tasks completed, before the call returns; thus significant time
!! can appear to be spent in this subroutine.
!!
!! @note
!! An actual, non-stub implementation of this interface is only available if
!! Flash-X is configured and linked appropriately.
!! See Orchestration_setupPipelineForGpuTasks for more information.
!!
!! @param nThreads              Number of threads in team to activate.
!!                              Used for checking that the configuration of
!!                              the orchestration runtime matches expectations.
!! @param nTilesPerPacket       The maximum number of tiles allowed in each
!!                              packet.
!!                              Used for checking that the configuration of
!!                              the orchestration runtime matches expectations.
!!                              Note that the presense of this non-optional
!!                              argument distinguishes the specific interface
!!                              Orchestration_teardownPipelineForGpuTasks from
!!                              Orchestration_teardownPipelineForCpuTasks.
subroutine Orchestration_teardownPipelineForGpuTasks(nThreads, nTilesPerPacket)
   implicit none

   integer,                                intent(IN) :: nThreads
   integer,                                intent(IN) :: nTilesPerPacket

end subroutine Orchestration_teardownPipelineForGpuTasks

