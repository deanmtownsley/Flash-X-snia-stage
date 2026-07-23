!> @copyright Copyright 2024 UChicago Argonne, LLC and contributors
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
!> @ingroup TimeAdvance
!!
!! @brief Initialize any variables that are needed by TimeAdvance unit
!!
!! @details
!!
!! This procedure implements initializations that the TimeAdvance
!! unit needs.

subroutine TimeAdvance_init()
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use TimeAdvance_data, ONLY: ta_useAsyncGrav
  implicit none

  call RuntimeParameters_get("ta_useAsyncGrav", ta_useAsyncGrav)
end subroutine TimeAdvance_init
