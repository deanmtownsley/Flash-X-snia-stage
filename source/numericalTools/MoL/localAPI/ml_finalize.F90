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
!! @brief ml_finalize stub

!> @ingroup MoLPrivate
!!
!! @brief Finalize a MoL implementation unit
!!
!! @details
!! @anchor ml_finalize_stub
!!
!! This procedure is responsible deallocating all memory allocated
!! by a specific implementation.  This procedure will be called once
!! during `MoL_finalize`, and must be matched by a call to `ml_init`
!! during `MoL_init`
subroutine ml_finalize()
   implicit none

   return
end subroutine ml_finalize
