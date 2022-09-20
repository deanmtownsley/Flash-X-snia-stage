!> @file source/numericalTools/MoL/MoL_regrid.F90
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
!! @brief MoL_regrid stub
!! @ingroup MoL

!> @brief Regrid MoL data structures
!! @anchor MoL_regrid_stub
!!
!! @pre MoL has been initialized
!! @pre Grid has been initialized
!! @pre Evolved variables are registered with MoL
!!
!! @post The dimensions and refinement of MoL's data structures will match
!!       those of Grid
!!
!! @sideeffect MoL's data structures may be reallocated
!!
!! @returns None
!!
!! @ingroup MoL
subroutine MoL_regrid()
   implicit none

   return
end subroutine MoL_regrid
