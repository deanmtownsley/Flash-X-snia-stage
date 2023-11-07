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
!! @brief Spacetime_unitTest stub

!> @ingroup physics_Spacetime
!!
!! @brief Perform a unit test of a Spacetime implementation
!!
!! @details
!! @anchor Spacetime_unitTest_stub
!!
!! @param fileUnit  The file unit to write the results to
!! @param perfect   True if no errors occurred; false if errors occurred
subroutine Spacetime_unitTest(fileUnit, perfect)
   implicit none

   integer, intent(in) :: fileUnit
   logical, intent(inout) :: perfect

   perfect = .false.

   return
end subroutine Spacetime_unitTest
