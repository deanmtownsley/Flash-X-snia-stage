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
!! @brief Spacetime_molPostRegrid stub

!> @ingroup physics_Spacetime
!!
!! @brief Perform any necessary work after any re-gridding operation
!!
!! @details
!! @anchor Spacetime_molPostRegrid_stub
!!
!! This procedure is responsible for performing any required tasks
!! after any changes to the evolved variables during a re-gridding
!! operation.  Some tracked quantities may need to be recalculated
!! from the new values of the evolved variables.
!!
!! @param  t  Time of the current solution of the evolved variables
subroutine Spacetime_molPostRegrid(t)
   implicit none

   real, intent(in) :: t

   return
end subroutine Spacetime_molPostRegrid
