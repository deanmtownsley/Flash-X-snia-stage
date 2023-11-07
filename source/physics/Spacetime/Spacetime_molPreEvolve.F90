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
!! @brief Spacetime_molPreEvolve stub

!> @ingroup physics_Spacetime
!!
!! @brief Perform any necessary work prior to starting evolution
!!
!! @details
!! @anchor Spacetime_molPreEvolve_stub
!!
!! This procedure is responsible for performing any required tasks
!! prior to the start of the evolution time loop.  Potential work
!! performed here will including setting/updating evolved variables
!! based on the initial data provided during `Simulation_initBlock`,
!! and anything else that requires all other units to have completed
!! initialization
!!
!! @param  t  Time at the start of the evolution
subroutine Spacetime_molPreEvolve(t)
   implicit none

   real, intent(in) :: t

   return
end subroutine Spacetime_molPreEvolve
