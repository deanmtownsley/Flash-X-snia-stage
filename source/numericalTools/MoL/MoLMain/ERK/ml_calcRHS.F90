!!****if* source/numericalTools/MoL/MoLMain/ERK/ml_calcRHS
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
subroutine ml_calcRHS(rhsType, rhsStruct, t)
   use MoL_functions, only: MoL_rhsE, MoL_rhsI, MoL_rhsF

   use ml_memInterface, only: ml_memSetActiveRHS, ml_memReleaseActiveRHS, ml_memZero

   implicit none

   integer, intent(in) :: rhsType, rhsStruct
   real, intent(in) :: t

   ! Zero-out RHS memory
   call ml_memZero(rhsStruct)

   call ml_memSetActiveRHS(rhsStruct)

   call MoL_rhsE(t)
   call MoL_rhsI(t)
   call MoL_rhsF(t)

   call ml_memReleaseActiveRHS()
end subroutine ml_calcRHS
