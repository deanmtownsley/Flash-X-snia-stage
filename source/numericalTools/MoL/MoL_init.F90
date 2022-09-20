!> @file source/numericalTools/MoL/MoL_init.F90
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
!! @brief MoL_init stub
!! @ingroup MoL

!> @brief Initialize the method of lines unit
!! @anchor MoL_init_stub
!!
!! @pre The MoL unit has not been initialized
!! @post All MoL runtime parameters have been read and saved, and the unit
!!       is ready to receive evolved variable registrations
!!
!! @sideeffect MoL unit-wide evolved-variable tracking and integeration scheme
!!             internal variables have been allocated and initialized
!!
!! @returns None
!!
!! @ingroup MoL
subroutine MoL_init()
   implicit none

   return
end subroutine MoL_init
