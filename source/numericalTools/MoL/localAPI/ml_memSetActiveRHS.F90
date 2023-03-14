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
!! @brief ml_memSetActiveRHS stub

!> @ingroup MoLPrivate
!!
!! @brief Set the active RHS data structure
!!
!! @details
!! @anchor ml_memSetActiveRHS_stub
!!
!! This procedure will set the active RHS term that is associated with
!! requests for the MOL_RHS data pointer.  All subsequent requests for
!! pointers to MOL_RHS will now target the pointer torwards the active RHS
!!
!! @param irhs Identifier of the RHS data structure to associate with MOL_RHS
subroutine ml_memSetActiveRHS(irhs)

   implicit none

   integer, intent(in) :: irhs

   return
end subroutine ml_memSetActiveRHS
