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
!! @brief Spacetime_finalize stub

!> @ingroup physics_Spacetime
!!
!! @brief Finalize the Spacetime unit
!!
!! @details
!! @anchor Spacetime_finalize_stub
!!
!! This procedure must be called once and only once during shutdown, and
!! must be matched by a call to @ref spacetime_init during startup.  Any
!! memory allocated by this unit will be freed.
subroutine Spacetime_finalize()
   implicit none

   return
end subroutine Spacetime_finalize
