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

!> @ingroup OrchestrationMilhoja
!! @stubref{Orchestration_finalize}
!!
!! @brief Concrete implementation of Orchestration_finalize
subroutine Orchestration_finalize()
    use milhoja_types_mod,       ONLY : MILHOJA_INT
    use milhoja_runtime_mod,     ONLY : milhoja_runtime_finalize

    use Orchestration_interface, ONLY : Orchestration_checkInternalError

    implicit none

    integer(MILHOJA_INT) :: MH_ierr

    CALL milhoja_runtime_finalize(MH_ierr)
    CALL Orchestration_checkInternalError("Orchestration_finalize", MH_ierr)
end subroutine Orchestration_finalize

