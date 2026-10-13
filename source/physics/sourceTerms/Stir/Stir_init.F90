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
!! @brief Stir_init stub

!> @ingroup Stir
!!
!! @brief Initialize turbulence stirring module
!!
!! @details
!! @anchor Stir_init_stub
!!
!! Does initializations that the Stir unit might need.
!!
!! @param restart  boolean flag to indicate whether restarting from checkpoint file

subroutine Stir_init(restart)
  implicit none
  logical, intent(in) :: restart
  return
end subroutine Stir_init
