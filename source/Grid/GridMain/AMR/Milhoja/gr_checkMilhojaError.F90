#include "Milhoja_interface_error_codes.h"

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
!! A helper routine that checks if the Milhoja call that returned the given
!! error code ended successfully or in failure.  If failure, then this routine
!! composes an informative error message, which it passes to Driver_abort.
!!
!! @param routineName  The name of the Flash-X routine that called this routine.
!!                     This information is included in the error message.
!! @param MH_errorCode The Milhoja error code to check.
subroutine gr_checkMilhojaError(routineName, MH_errorCode)
    use milhoja_types_mod,  ONLY : MILHOJA_INT
    use milhoja_errors_mod, ONLY : MAX_ERROR_LENGTH, &
                                   milhoja_errorMessage

    use Driver_interface, ONLY : Driver_abort

    implicit none

    character(LEN=*),     intent(IN) :: routineName
    integer(MILHOJA_INT), intent(IN) :: MH_errorCode

    character(LEN=  MAX_ERROR_LENGTH) :: msg
    character(LEN=3*MAX_ERROR_LENGTH) :: buffer

    if (MH_errorCode /= MILHOJA_SUCCESS) then
        CALL milhoja_errorMessage(MH_errorCode, msg)

        write(buffer,'(3A,I0,2A)') "[", TRIM(routineName), "] Milhoja error ", &
                                   MH_errorCode, " - ", TRIM(msg)

        CALL Driver_abort(buffer)
    end if
end subroutine gr_checkMilhojaError

