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
!! @brief MoL_advance stub

!> @ingroup MoL
!!
!! @brief Take a timestep from t to t+dt
!!
!! @details
!! @anchor MoL_advance_stub
!!
!! This procedure is the main entry point into the MoL integrator,
!! and should be called once per time step.  All evolved variables
!! will be advanced from the time `t` to the time `t + dt` based
!! on the evaluation of the RHS terms for the variables' equations.
!!
!! @todo When/if subcyling is available in Flash-X, this will extend
!!       to include a level-indicator as well
!!
!! @param t  Current time
!! @param dt Size of the timestep to take
subroutine MoL_advance(t, dt)
   implicit none

   real, intent(in) :: t, dt

   return
end subroutine MoL_advance
