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
!! @brief MoL_getRHSIndex implementation

!> @ingroup MoLMain
!!
!! @brief Implements MoL_getRHSIndex
!!
!! @stubref{MoL_getRHSIndex}
function MoL_getRHSIndex(evolIndex) result(rhsIndex)
   use ml_variables, only: ml_unk_to_scratch
   implicit none

   integer :: rhsIndex

   integer, intent(in) :: evolIndex

   rhsIndex = ml_unk_to_scratch(evolIndex)
end function MoL_getRHSIndex
