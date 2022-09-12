!!****if* source/numericalTools/MoL/MoLMain/MoL_registerRHS
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
!!      MoL_registerRHS
!!
!!  SYNOPSIS
!!
!!      call MoL_registerRHS(integer,   intent(in) :: rhsType,
!!                           procedure,            :: rhsFunc)
!!
!!  DESCRIPTION
!!
!!      Register a function with MoL to evaluate one of the three
!!      right-hand side types
!!
!! ARGUMENTS
!!
!!      rhsType : One of the following function types defined in MoL.h
!!                   MOL_RHS_EXPLICIT  -  RHS for (slow) explicit terms
!!                   MOL_RHS_IMPLICIT  -  RHS for (slow) implicit terms
!!                   MOL_RHS_FAST      -  RHS for (fast) explicit terms
!!      rhsFunc : Procedure to register
!!
!!***
subroutine MoL_registerRHS(rhsType, rhsFunc)
   use ml_functions, only: ml_rhs_t, ml_rhsE, ml_rhsI, ml_rhsF

   use ml_interface, only: ml_error

#include "MoL.h"

   implicit none

   integer, intent(in) :: rhsType
   procedure(ml_rhs_t) :: rhsFunc

   select case (rhsType)
   case (MOL_RHS_EXPLICIT)
      if (associated(ml_rhsE)) nullify (ml_rhsE)
      ml_rhsE => rhsFunc

   case (MOL_RHS_IMPLICIT)
      if (associated(ml_rhsI)) nullify (ml_rhsI)
      ml_rhsI => rhsFunc

   case (MOL_RHS_FAST)
      if (associated(ml_rhsF)) nullify (ml_rhsF)
      ml_rhsF => rhsFunc

   case default
      call ml_error("Attempting to register unknown RHS function type")
   end select ! rhsType
end subroutine MoL_registerRHS
