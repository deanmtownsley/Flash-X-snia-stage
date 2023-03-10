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
!! @brief MoL_finalize stub

!> @ingroup MoL
!!
!! @brief Finalize the method of lines unit
!!
!! @details
!! @anchor MoL_finalize_stub
!!
!! This procedure will deallocate memory managed by MoL and call
!! all implementation specific finalization procedures.  A single call
!! to this procedure during program execution finalization is required
!! if a call to `MoL_initialize` has been made
subroutine MoL_finalize()
   implicit none

   return
end subroutine MoL_finalize
