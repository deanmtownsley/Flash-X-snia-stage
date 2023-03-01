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
!! @brief MoL_getRHSIndex stub

!> @ingroup MoL
!!
!! @brief Obtain the RHS index for a registered evolved variable
!!
!! @details
!! @anchor MoL_getRHSIndex_stub
!!
!! This function will provide the RHS index of a registered evolved variable
!! as determined during calls to @ref MoL_registerVariable.
!!
!! @note This RHS index is also provided as an output of @ref MoL_registerVariable
!!
!! @param evolIndex  Index of the evolved variable in UNK
!! @param rhsIndex   Output index of the evolved variable in MoL RHS
function MoL_getRHSIndex(evolIndex) result(rhsIndex)
   implicit none

   integer :: rhsIndex

   integer, intent(in) :: evolIndex

   rhsIndex = -1

   return
end function MoL_getRHSIndex
