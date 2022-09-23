!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! @brief ml_memReleaseActiveRHS stub

!> @ingroup MoLPrivate
!!
!! @brief Releasse the active RHS data structure
!!
!! @details
!! @anchor ml_memReleaseActiveRHS_stub
!!
!! This procedure will unset the active RHS term that is associated with
!! requests for the MOL_RHS data pointer
!!
!! @post MOL_RHS will no longer refer to a specific stage/type of RHS
!!       storage in MoL's scratch memory
!!
!! @returns None
subroutine ml_memReleaseActiveRHS(irhs)

   implicit none

   integer, intent(in) :: irhs

   return
end subroutine ml_memReleaseActiveRHS
