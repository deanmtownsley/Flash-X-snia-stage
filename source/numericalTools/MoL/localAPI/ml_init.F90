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
!! @brief ml_init stub

!> @ingroup MoLPrivate
!!
!! @brief Initialize a MoL implementation unit
!!
!! @details
!! @anchor ml_init_stub
!!
!! This procedure is responsible for reading and runtime parameters,
!! allocating memory, and all other tasks that are specific to the
!! implementation.  It will be called only once during `MoL_init`, and
!! will be matched by a call to `ml_finalize` during `MoL_finalize`
subroutine ml_init()
   implicit none

   return
end subroutine ml_init
