!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! @brief ml_calcRHS implementation for FBE

!> @ingroup MoLFBE
!!
!! @brief Implements ml_calcRHS for FBE
!!
!! @stubref{ml_calcRHS}
subroutine ml_calcRHS(rhsType, rhsStruct, t)
   use ml_functions, only: ml_rhsE, ml_rhsF

   use ml_memInterface, only: ml_memZero

   implicit none

   integer, intent(in) :: rhsType, rhsStruct
   real, intent(in) :: t

   ! Zero-out RHS memory
   call ml_memZero(rhsStruct)

   ! No need to set an active RHS - only MOL_RHS exists
   ! Both explicit and fast-explicit will be added here
   call ml_rhsE(t)
   call ml_rhsF(t)

end subroutine ml_calcRHS
