!> @file
!!
!! @copyright Copyright 2023 UChicago Argonne, LLC and contributors
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
!! @brief MoL_releaseFunctions stub

!> @ingroup MoL
!!
!! @brief Release all registered functions
!!
!! @details
!! @anchor MoL_releaseFunctions_stub
!!
!! Releasing registered functions will reset MoL's internal procedure pointers
!! to target their default "stub" implementations
subroutine MoL_releaseFunctions()
   implicit none

   return
end subroutine MoL_releaseFunctions
