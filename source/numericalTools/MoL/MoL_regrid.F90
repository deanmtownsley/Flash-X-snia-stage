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
!! @brief MoL_regrid stub

!> @ingroup MoL
!!
!! @brief Regrid MoL data structures
!!
!! @details
!! @anchor MoL_regrid_stub
!!
!! This procedure will update the shape of MoL's data structures to
!! match the current grid.  This may require (re)allocation of the
!! data structures in some implementations.  Space will be allocated
!! for all evolved variables in every data structure.
subroutine MoL_regrid()
   implicit none

   return
end subroutine MoL_regrid
