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
!! @stubref{Orchestration_pushTileToPipeline}
!!
!! @brief Stub implementation of Orchestration_pushTileToPipeline

subroutine Orchestration_pushTileToPipeline(prototype_Cptr, nThreads, &
                                            tileCInfo)
    use iso_c_binding, ONLY : C_PTR, c_loc

    use Orchestration_interfaceTypeDecl, ONLY: Orchestration_tileCInfo_t

    implicit none

    type(C_PTR),                            intent(IN) :: prototype_Cptr
    integer,                                intent(IN) :: nThreads
    type(Orchestration_tileCInfo_t),intent(IN),target :: tileCInfo

end subroutine Orchestration_pushTileToPipeline
! Local Variables:
! f90-program-indent: 4
! f90-do-indent: 3
! f90-type-indent: 3
! indent-tabs-mode: nil
! End:
