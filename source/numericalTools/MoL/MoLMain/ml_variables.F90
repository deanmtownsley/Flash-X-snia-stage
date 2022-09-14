!!****if* source/numericalTools/MoL/MoLMain/ml_variables
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
!!      ml_variables
!!
!!  SYNOPSIS
!!
!!      use ml_variables, only:
!!
!!  DESCRIPTION
!!
!!      Evolved-variable tracking used internally by MoL
!!
!!***

module ml_variables

   implicit none

   type :: ml_variable_t
      character(len=:), allocatable :: name
      integer :: evolIndex, rhsIndex
   end type ml_variable_t

   type(ml_variable_t), allocatable, save :: ml_vars(:)
   integer, save :: ml_nvars = 0

   ! For convenience
   integer, allocatable, save :: ml_unk_mask(:), ml_scratch_mask(:)

end module ml_variables
