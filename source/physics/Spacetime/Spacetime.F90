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
!! @brief Spacetime stub

!> @ingroup physics_Spacetime
!!
!! @brief Evolve the spacetime variables by the provided timestep
!!
!! @details
!! @anchor Spacetime_stub
!!
!! This procedure is intended to evolve the spacetime variables when not
!! utilizing the @ref MoL "MoL Unit".  It is not advisable to use this
!! procedure along side other physics units that are evolving equations
!! that may tightly couple to those evolved by the Spacetime unit.  It
!! is not guaranteed that an implementation of Spacetime will provide a
!! non-stub implementation of this procedure.
!!
!! @param  t   The time at the start of the timestep
!! @param  dt  The size of the timestep
subroutine Spacetime(t, dt)
   implicit none

   real, intent(in) :: t, dt

   return
end subroutine Spacetime
