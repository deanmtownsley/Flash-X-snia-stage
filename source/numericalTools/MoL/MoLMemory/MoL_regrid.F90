!> @file source/numericalTools/MoL/MoLMemory/MoL_regrid.F90
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
!! @brief MoL_regrid implementation
!! @ingroup MoLMemory

!> @brief Implements MoL_regrid
!! @details This implementation will work will all Grid backends
!!
!! @ref MoL_regrid_stub "See stub documentation"
!!
!! @ingroup MoLMemory
subroutine MoL_regrid()
   use ml_memInterface, only: ml_memAlloc, ml_memFree

   implicit none

   logical, save :: first = .true.

   ! Only need to do this once for UG/Paramesh
   if (.not. first) return

   call ml_memFree
   call ml_memAlloc

   first = .false.
end subroutine MoL_regrid
