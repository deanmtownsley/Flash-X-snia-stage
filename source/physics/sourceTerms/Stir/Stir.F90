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
!! @brief Stir stub

!> @ingroup Stir
!!
!! @brief Apply turbulence stirring
!!
!! @details
!! @anchor Stir_stub
!!
!! Apply turbulence stirring via TurbGen (see Federrath et al. 2010; 2022).
!!
!! @param dt    current time step
!! @param pass  optional in case of split time advancement to indicate 1st or 2nd half step

subroutine Stir(dt, pass)
  implicit none
  real,intent(IN) :: dt
  integer, intent(in), optional :: pass
  return
end subroutine Stir
