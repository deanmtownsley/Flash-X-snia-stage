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
!! @brief ml_memFree stub

!> @ingroup MoLPrivate
!!
!! @brief Dellocate MoL scratch memory data structures
!!
!! @details
!! @anchor ml_memFree_stub
!!
!! This procedure will be call immediately before allocation to ensure that
!! memory is not already allocated, and upon finalization of the unit
subroutine ml_memFree()
   implicit none

   return
end subroutine ml_memFree
