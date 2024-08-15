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
!! @stubref{Orchestration_pushTileToCpuGpuPipeline}
!!
!! @brief Concrete implementation of Orchestration_pushTileToCpuGpuPipeline
#include "Milhoja.h"
subroutine Orchestration_pushTileToCpuGpuPipeline(pktProto_Cptr, tileProto_Cptr, &
                                            nThreads, &
                                            tileCInfo)
    use iso_c_binding, ONLY : C_PTR, c_loc

    use milhoja_types_mod,   ONLY : MILHOJA_INT
#if ! defined(RUNTIME_MUST_USE_TILEITER) && defined(RUNTIME_SUPPORT_DATAPACKETS)
    use milhoja_runtime_mod, ONLY : milhoja_runtime_pushTileToCpuGpuPipeline
#endif

    use Driver_interface,        ONLY : Driver_abort
    use Orchestration_interfaceTypeDecl, ONLY: Orchestration_tileCInfo_t
    use Orchestration_interface, ONLY : Orchestration_checkInternalError

    implicit none

    type(C_PTR),                            intent(IN) :: pktProto_Cptr
    type(C_PTR),                            intent(IN) :: tileProto_Cptr
    integer,                                intent(IN) :: nThreads
    type(Orchestration_tileCInfo_t),intent(IN),target :: tileCInfo

    integer(MILHOJA_INT) :: MH_nThreads
    integer(MILHOJA_INT) :: MH_ierr
    type(C_PTR) :: MH_tileCInfo_Cp

    MH_nThreads = INT(nThreads, kind=MILHOJA_INT)
    MH_tileCInfo_Cp = c_loc(tileCInfo)

#if ! defined(RUNTIME_MUST_USE_TILEITER) && defined(RUNTIME_SUPPORT_DATAPACKETS)
    CALL milhoja_runtime_pushTileToCpuGpuPipeline(pktProto_Cptr, &
                                                  tileProto_Cptr, &
                                          MH_nThreads, MH_tileCInfo_Cp, MH_ierr)
    CALL Orchestration_checkInternalError("Orchestration_pushTileToCpuGpuPipeline", MH_ierr)
#else
    CALL Driver_abort("Orchestration_pushTileToCpuGpuPipeline: milhoja_runtime_pushTileToCpuGpuPipeline disabled")
#endif
end subroutine Orchestration_pushTileToCpuGpuPipeline
! Local Variables:
! f90-program-indent: 4
! f90-do-indent: 3
! f90-type-indent: 3
! indent-tabs-mode: nil
! End:
