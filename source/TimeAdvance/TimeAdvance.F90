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
!> @ingroup TimeAdvance
!!
!! @brief Take a timestep from t-dt to t
!!
!! @details
!! @anchor TimeAdvance_stub
!!
!! This procedure is the main entry point into the  integrator,
!! and should be called once per time step.  All evolved variables
!! will be advanced from the time `t-dt` to the time `t` 
!!
!!
!! @param time  Current time
!! @param dt Size of the timestep to take
!! @param dtold Size of the previous timestep 

subroutine TimeAdvance(dt, dtold, time)
  real, intent(IN) :: dt, dtold, time
end subroutine TimeAdvance
