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
!! @stubref{Orchestration_pushTileToExtGpuPipeline}
!!
!! @brief Stub implementation of Orchestration_pushTileToExtGpuPipeline
!!
!! @details
!! Push a tile to the orchestration runtime for processing by the task
!! functions that were passed when the runtime system was set up with a
!! preceding Orchestration_setupPipelineForExtGpuTasks call.
!! This variant of the push interface is for the Extended GPU thread team configuration.
!! A sequence of Orchestration_pushTileToExtGpuPipeline calls (usually
!! occurring in an iterator loop over tiles) has to be bracketed by calls for
!! setting up and for tearing down the desired thread team configuration with
!! the desired task functions, such that the sequence of Orchestration calls is
!! - `call Orchestration_setupPipelineForExtGpuTasks`
!! - `call Orchestration_pushTileToExtGpuPipeline` [...]
!! - `call Orchestration_teardownPipelineForExtGpuTasks`
!!
!! @param prototype_Cptr        Pointer to a prototype datapacket, from which
!!                              the orchestration runtime will generate new
!!                              datapackets.
!!                              Prototype datapackets are created by calling
!!                              a function like, for example,
!!                              instantiate_gpu_tf_hydro_packet_C in the case
!!                              of orchestrating work for a Hydro operation.
!! @param nThreads              The number of threads to activate in the teams
!!                              that apply the task functions.
!!                              Used here only for checking that the current
!!                              configuration of the orchestration runtime,
!!                              established by the previous call of
!!                              Orchestration_setupPipelineForExtGpuTasks,
!!                              matches expectations.
!! tileCInfo                    An object of C-compatible derived type holding
!!                              information that identifies and describes the
!!                              tile on which work is to be done; including
!!                              pointers to the Grid unit's raw (real) data
!!                              associated with the tile.
subroutine Orchestration_pushTileToExtGpuPipeline(prototype_Cptr, nThreads, &
                                            tileCInfo)
    use iso_c_binding, ONLY : C_PTR
    use Orchestration_interfaceTypeDecl, ONLY: Orchestration_tileCInfo_t

    implicit none

    type(C_PTR),                            intent(IN) :: prototype_Cptr
    integer,                                intent(IN) :: nThreads
    type(Orchestration_tileCInfo_t),intent(IN),target :: tileCInfo

end subroutine Orchestration_pushTileToExtGpuPipeline
! Local Variables:
! f90-program-indent: 4
! f90-do-indent: 3
! f90-type-indent: 3
! indent-tabs-mode: nil
! End:
