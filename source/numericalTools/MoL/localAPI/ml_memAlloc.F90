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
!! @brief ml_memAlloc stub

!> @ingroup MoLPrivate
!!
!! @brief Allocate MoL scratch memory data structures
!!
!! @details
!! @anchor ml_memAlloc_stub
!!
!! This procedure may be used either at startup or after regridding.  Space
!! for all evolved variables will be allocated
subroutine ml_memAlloc()
   implicit none

   return
end subroutine ml_memAlloc
