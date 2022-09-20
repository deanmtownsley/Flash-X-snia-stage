!> @file source/numericalTools/MoL/MoLMain/MoL_finalize.F90
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
!! @brief MoL_finalize implementation
!! @ingroup MoLMain

!> @brief Implements MoL_finalize
!!
!! @ref MoL_finalize_stub "See stub documentation"
!!
!! @ingroup MoLMain
subroutine MoL_finalize()
   use ml_interface, only: ml_finalize
   use ml_memInterface, only: ml_memFree

   implicit none

   call ml_memFree
   call ml_finalize
end subroutine MoL_finalize
