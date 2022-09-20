!> @file source/numericalTools/MoL/MoL_releaseFunctions.F90
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
!! @brief MoL_releaseFunctions stub
!! @ingroup MoL

!> @brief Release all registered functions
!! @anchor MoL_releaseFunctions_stub
!!
!! @pre  MoL's procedure pointers target either registered or default procedures
!! @post MoL's procedure pointers will target default versions
!!
!! @sideeffect The target of MoL's procedure pointers will change
!!
!! @returns None
!!
!! @ingroup MoL
subroutine MoL_releaseFunctions()
   implicit none

   return
end subroutine MoL_releaseFunctions
