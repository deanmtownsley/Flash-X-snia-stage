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
!! @brief MoL_init stub

!> @ingroup MoL
!!
!! @brief Initialize the method of lines unit
!!
!! @details
!! @anchor MoL_init_stub
!!
!! This procedure should not be called more than once per program execution
!! during startup, only after RuntimeParameters has been initialized, but before
!! Simulation and all physics units are initiailized.  All of MoL's runtime parameters
!! will be processed, and the unit will be prepared to receive variable registrations.
!! All implementation-specific intialization procedures will be called subequently during
!! execution of this procedure.
subroutine MoL_init()
   implicit none

   return
end subroutine MoL_init
