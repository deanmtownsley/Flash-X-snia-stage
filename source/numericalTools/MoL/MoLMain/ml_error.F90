!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief ml_error implementation

!> @ingroup MoLMain
!!
!! @brief Implements ml_error
!!
!! @stubref{ml_error}
subroutine ml_error(msg)
   use Driver_interface, only: Driver_abort

#include "constants.h"

   implicit none

   character(len=*), intent(in) :: msg

   ! Error-messaging is always turned on
   call Driver_abort("[MoL] ERROR: "//msg)
end subroutine ml_error
