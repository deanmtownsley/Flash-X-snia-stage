!!****if* source/numericalTools/MoL/MoLMain/MR/ml_calcRHS
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!!  NAME
!!
!!      ml_calcRHS
!!
!!  SYNOPSIS
!!
!!      call ml_calcRHS(integer, intent(in) :: rhsType
!!                      integer, intent(in) :: rhsStruct
!!                      real,    intent(in) :: t)
!!
!!  DESCRIPTION
!!
!!      Calculate specified RHS type and store in  requested struct
!!
!!  ARGUMENTS
!!
!!      rhsType   : The type of RHS, one of:
!!                  - MOL_RHS_EXPLICIT
!!                  - MOL_RHS_IMPLICIT
!!                  - MOL_RHS_FAST
!!      rhsStruct : MoL memory data-struct to store RHS in
!!      t         : The time of the RHS is to be evaluated at
!!
!!***
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
!! @brief ml_calcRHS implementation for MR

!> @ingroup MoLMR
!!
!! @brief Implements ml_calcRHS for MR
!!
!! @stubref{ml_calcRHS}
subroutine ml_calcRHS(rhsType, rhsStruct, t)
   use ml_functions, only: ml_rhsE, ml_rhsI, ml_rhsF
   use ml_memInterface, only: ml_memSetActiveRHS, ml_memReleaseActiveRHS, ml_memZero

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

   implicit none

   integer, intent(in) :: rhsType, rhsStruct
   real, intent(in) :: t

   ! Zero-out RHS memory
   call ml_memZero(rhsStruct)

   call ml_memSetActiveRHS(rhsStruct)

   select case (rhsType)
   case (MOL_RHS_EXPLICIT)
      call ml_rhsE(t)

   case (MOL_RHS_IMPLICIT)
      call ml_rhsI(t)

   case (MOL_RHS_FAST)
      call ml_rhsF(t)
   end select

   call ml_memReleaseActiveRHS()

end subroutine ml_calcRHS
