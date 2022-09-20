!> @file source/numericalTools/MoL/MoL_finalize.F90
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
!! @brief MoL_finalize stub
!! @ingroup MoL

!> @brief Finalize the method of lines unit
!! @anchor MoL_finalize_stub
!!
!! @pre The MoL has been initialized
!! @post Evolved variables are no longer registered and MoL memory structures
!!       are inaccessible
!!
!! @sideeffect All MoL allocated memory will be freed
!!
!! @returns None
!!
!! @ingroup MoL
subroutine MoL_finalize()
   implicit none

   return
end subroutine MoL_finalize
