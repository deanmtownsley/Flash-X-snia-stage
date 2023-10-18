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
!! @brief Spacetime_init stub

!> @ingroup physics_Spacetime
!!
!! @brief Initialize the Spacetime unit
!!
!! @details
!! @anchor Spacetime_init_stub
!!
!! This procedure must be called once and only once during startup, and must
!! be called only after the @ref MoL "MoL" unit has been initialized.  Tasks
!! performed by this procedure will include reading runtime parameters,
!! registering evolved variables, and any other initialization tasks required
!! by the chosen implementation. A call to this procedure must be matched by
!! a call to @ref spacetime_finalize
subroutine Spacetime_init()
   implicit none

   return
end subroutine Spacetime_init
