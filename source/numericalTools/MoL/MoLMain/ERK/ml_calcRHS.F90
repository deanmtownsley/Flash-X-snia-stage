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
!! @brief ml_calcRHS implementation for ERK

!> @ingroup MoLERK
!!
!! @brief Implements ml_calcRHS for ERK
!!
!! @stubref{ml_calcRHS}
subroutine ml_calcRHS(rhsType, rhsStruct, t, dtWeight)
   use ml_functions, only: ml_rhsE, ml_rhsI, ml_rhsF

   use ml_memInterface, only: ml_memZero

   implicit none

   integer, intent(in) :: rhsType, rhsStruct
   real, intent(in) :: t
   real, intent(in) :: dtWeight

   ! Zero-out RHS memory
   call ml_memZero(rhsStruct)

   call ml_rhsE(t, rhsStruct, dtWeight)
   call ml_rhsI(t, rhsStruct, dtWeight)
   call ml_rhsF(t, rhsStruct, dtWeight)
end subroutine ml_calcRHS
