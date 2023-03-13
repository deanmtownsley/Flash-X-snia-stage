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
!!
!! @brief ml_memZero stub

!> @ingroup MoLPrivate
!!
!! @brief Initialize a MoL data structure to zero
!!
!! @details
!! @anchor ml_memZero_stub
!!
!! @param dataStruct  Destintation data structure
subroutine ml_memZero(dataStruct)
   implicit none

   integer, intent(in) :: dataStruct

   return
end subroutine ml_memZero
