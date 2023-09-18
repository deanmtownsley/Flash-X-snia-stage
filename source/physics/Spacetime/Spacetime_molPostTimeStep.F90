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
!! @brief Spacetime_molPostTimeStep stub

!> @ingroup physics_Spacetime
!!
!! @brief Perform any necessary work after a full timestep
!!
!! @details
!! @anchor Spacetime_molPostTimeStep_stub
!!
!! This procedure is responsible for performing any required tasks
!! after the evolved variables are updated by a full timestep.  Some of
!! these tasks may include analytic and diagnostic calculations (e.g. the
!! ADM mass)
!!
!! @param  t  Time of the current solution of the evolved variables
subroutine Spacetime_molPostTimeStep(t)
   implicit none

   real, intent(in) :: t

   return
end subroutine Spacetime_molPostTimeStep
