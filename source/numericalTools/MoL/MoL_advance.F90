!> @file source/numericalTools/MoL/MoL_advance.F90
!!
!! @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! @brief MoL_advance stub
!! @ingroup MoL

!> @brief Take a timestep from t to t+dt
!! @anchor MoL_advance_stub
!!
!! @param t  Current time
!! @param dt Size of the timestep to take
!!
!! @pre The solution variables are at a time `t`
!! @post The solution variables will be advanced to a time `t+dt`
!!
!! @returns None
!!
!! @todo When/if subcyling is available in Flash-X, this will extend
!!       to include a level-indicator as well
!!
!! @ingroup MoL
subroutine MoL_advance(t, dt)
   implicit none

   real, intent(in) :: t, dt

   return
end subroutine MoL_advance
