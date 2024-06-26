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
!! @stubref{Orchestration_pushTileToGpuPipeline}
!!
!! @brief Stub implementation of Orchestration_pushTileToGpuPipeline

!! @param prototype_Cptr        Pointer to a prototype datapacket, from which
!!                              the orchestration runtime will generate new
!!                              datapackets.
!!                              Prototype datapackets are created by calling
!!                              a function like, for example,
!!                              instantiate_gpu_tf_hydro_packet_C in the case
!!                              of orchestrating work for a Hydro operation.
!! @param nThreads              The number of threads to activate in the team
!!                              that applies the task function.
!! tileCInfo                    An object of C-compatible derived type holding
!!                              information that identifies and describes the
!!                              tile on which work is to be done; including
!!                              pointers to the Grid unit's raw (real) data
!!                              associated with the tile.
subroutine Orchestration_pushTileToGpuPipeline(prototype_Cptr, nThreads, &
                                            tileCInfo)
   use iso_c_binding, ONLY : C_PTR

   use Orchestration_interfaceTypeDecl, ONLY: Orchestration_tileCInfo_t

   implicit none

   type(C_PTR),                            intent(IN) :: prototype_Cptr
   integer,                                intent(IN) :: nThreads
   type(Orchestration_tileCInfo_t),intent(IN),target :: tileCInfo
end subroutine Orchestration_pushTileToGpuPipeline
! Local Variables:
! f90-program-indent: 3
! f90-do-indent: 3
! f90-type-indent: 3
! indent-tabs-mode: nil
! End:
