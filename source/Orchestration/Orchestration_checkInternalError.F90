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

!> @ingroup Orchestration
!! @anchor Orchestration_checkInternalError_stub
!!
!! @brief Print an error message and abort if internal error reported
!!
!! @details
!! A helper routine that checks if the low-level Orchestration implementation
!! (e.g., Milhoja) call that returned the given error code ended successfully or in
!! failure.  If failure, then this routine composes an informative error message,
!! which it passes to Driver_abort.
!!
!! Since this routine deals with internal errors, it should be used only internally
!! and, therefore, should not be in the public interface.  However, this routine is
!! in the public interface so that tools in the offline toolchain that write
!! Flash-X Fortran code can perform error checking on internal calls.
!!
!! @note
!! This subroutine is available in the Orchestration unit's public interface,
!! whether the Milhoja Orchestration implementation is included or not. However,
!! the stub implementation will always abort when the Milhoja Orchestration
!! implementation is included, since in that case the Milhoja implementation
!! should override the stub.
!!
!! @param routineName  The name of the Flash-X routine that called this routine.
!!                     This information is included in the error message.
!! @param MH_errorCode The Milhoja error code to check.
subroutine Orchestration_checkInternalError(routineName, MH_errorCode)
    use Orchestration_interfaceTypeDecl, ONLY : MILHOJA_INT
    use Driver_interface,  ONLY : Driver_abort

    implicit none

    character(LEN=*),     intent(IN) :: routineName
    integer(MILHOJA_INT), intent(IN) :: MH_errorCode

#ifdef FLASHX_ORCHESTRATION_MILHOJA
    ! A lack of error checking should not be a quiet no-op.  Prefer abort to
    ! warning in the name of defensive programming.
    CALL Driver_abort("[Orchestration_checkInternalError] No error checking")
#endif
end subroutine Orchestration_checkInternalError

