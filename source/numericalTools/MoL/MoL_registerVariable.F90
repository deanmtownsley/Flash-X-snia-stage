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
!! @brief MoL_registerVariable stub

!> @ingroup MoL
!!
!! @brief Register an evolved variable with MoL
!!
!! @details
!! @anchor MoL_registerVariable_stub
!!
!! It is necessary to inform MoL of which variables in UNK are to be evolved
!! so that MoL knows how many RHS variables to allocate and which ones will
!! correspond to which evolved variables.
!!
!! @note Duplicate variable registrations will be ignored unless warnings are
!!       set to trigger runtime errors in MoL
!!
!! @todo This is intended as a temporary measure until a more suitable
!!       solution for MoL's scratch memory is decided.  The procedure name
!!       is intentionally similar to the same functionality provided
!!       in Cactus/ET as a reminder to myself find a better solution
!!
!! @param name       The name of the evolved variable
!! @param evolIndex  Index of the evolved variable in UNK
!! @param rhsIndex   Output index of the evolved variable in MoL RHS
!                    data structures
subroutine MoL_registerVariable(name, evolIndex, rhsIndex)
   implicit none

   character(len=*), intent(in) :: name
   integer, intent(in) :: evolIndex
   integer, intent(out) :: rhsIndex

   rhsIndex = -1

   return
end subroutine MoL_registerVariable
