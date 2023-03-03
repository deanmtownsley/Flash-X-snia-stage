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
!! @brief ml_status stub

!> @ingroup MoLPrivate
!!
!! @brief Issue a MoL-specific status message
!!
!! @details
!! @anchor ml_status_stub
!!
!! Status messages will only be displayed on the master process and if
!! MoL's verbosity level is set to `MOL_VERBOSITY_STATUS`
!!
!! @param msg  The status message to display
subroutine ml_status(msg)
   implicit none

   character(len=*), intent(in) :: msg

   return
end subroutine ml_status
