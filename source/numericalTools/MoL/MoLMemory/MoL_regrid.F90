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
!! @brief MoL_regrid implementation

!> @ingroup MoLMemory
!!
!! @brief Implements MoL_regrid
!!
!! @details
!! This implementation will work with UG and Paramesh.  The size of the
!! blocklists do not change after initial allocation, and the appropriate
!! block/tile refinement levels and layout will always be the same as
!! specified by the Grid iterators and tiles
!!
!! @stubref{MoL_regrid}
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
