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
!! @brief Spacetime_molPostFastUpdate stub

!> @ingroup physics_Spacetime
!!
!! @brief Perform any necessary work after an integration stage fast update
!!
!! @details
!! @anchor Spacetime_molPostFastUpdate_stub
!!
!! This procedure is responsible for performing any required tasks
!! after the evolved variables are updated during a fast integration stage
!!
!! @note When not using a multi-rate integrator, this procedure will be called
!!       after every update to the evolved variables.  Do no replicate
!!       the functionality provided here in or call this procedure from
!!       @ref spacetime_molpostupdate - doing so will end up doing the work
!!       performed here twice
!!
!! @param  t  Time of the current solution of the evolved variables
subroutine Spacetime_molPostFastUpdate(t)
   implicit none

   real, intent(in) :: t

   return
end subroutine Spacetime_molPostFastUpdate
