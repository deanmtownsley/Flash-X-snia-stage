!> @file source/numericalTools/MoL/MoL_registerVariable.F90
!!
!! @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! @brief MoL_registerVariable stub
!! @ingroup MoL

!> @brief Register an evolved variable with MoL
!! @anchor MoL_registerVariable_stub
!!
!! @param name       The name of the evolved variable
!! @param evolIndex  Index of the evolved variable in UNK
!! @param rhsIndex   Index of the evolved variable in MoL RHS data structures
!!
!! @pre `name` Is a unique identifier for the evolved variable
!! @pre `evolIndex` is valid variable index in UNK
!!
!! @post `rhsIndex` will be the index of the evolved variable in MoL's data structures
!! @post The number of registered variables will increase
!!
!! @sideeffect MoL internal variable trackers may be re-allocated
!!
!! @returns The index assigned to the evolved variable in MoL's data structures
!!
!! @details
!!    It is necessary to inform MoL of which variables in UNK are to be evolved
!!    so that MoL knows how many RHS variables to allocate and which ones will
!!    correspond to which evolved variables.
!!
!! @note Duplicate variable registrations will be ignored unless warnings are
!!       set to trigger runtime errors in MoL
!!
!!
!! @todo This is intended as a temporary measure until a more suitable
!!       solution for MoL's scratch memory is decided
!!
!! @ingroup MoL
subroutine MoL_registerVariable(name, evolIndex, rhsIndex)
   implicit none

   character(len=*), intent(in) :: name
   integer, intent(in) :: evolIndex
   integer, intent(out) :: rhsIndex

   rhsIndex = -1

   return
end subroutine MoL_registerVariable
