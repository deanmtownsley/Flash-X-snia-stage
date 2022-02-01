!!****f* source/Grid/Grid_addToVar
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
!! NAME
!!
!!  Grid_addToVar
!!
!! SYNOPSIS
!!
!!  call Grid_addToVar(integer(in) :: srcVar,
!!                     integer(in) :: destVar,
!!                     real(in) :: multFactor,
!!                     logical(in) :: reset)
!!
!! DESCRIPTION
!!   Compute solnData(srcVar,:,:,:)*multFactor and save in solnData(destVar,:,:,:).
!!
!!   If reset is true, the destination variable is first zeroed;
!!   otherwise the product is added to the existing values of destVar.
!!
!!   The operation is applied to interior cells of all LEAF blocks.
!!
!! ARGUMENTS
!!
!!
!!   srcVar : the state variables to be used in the RHS of the expression
!!
!!   destVar : the state variables to be used in the LHS of the expression
!!
!!   multFactor : multiplication factor
!!
!!   reset : indicates whether the destination variable should be zeroed first
!!
!! NOTES
!!
!!   srcVar == destVar is allowed and behaves as expected iff reset is .FALSE.
!!
!!   For a copy call Grid_addToVar(srcVar, destVar, 1.0, .true.)
!!
!!***

subroutine Grid_addToVar(srcVar, destVar, multFactor, reset)

  implicit none
  integer, intent(in) :: srcVar, destVar
  real,  intent(in) :: multFactor
  logical, intent(in) :: reset

end subroutine Grid_addToVar
