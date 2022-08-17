!!****if* source/numericalTools/MoL/MoLMain/MoL_variables
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
!!      MoL_variables
!!
!!  SYNOPSIS
!!
!!      use MoL_variables
!!
!!  DESCRIPTION
!!
!!      Evolved-variable tracking used internally by MoL
!!
!!***

module MoL_variables

   implicit none

   type :: MoL_variable_t
      character(len=:), allocatable :: name
      integer :: evolIndex, rhsIndex
   end type MoL_variable_t

   type(MoL_variable_t), allocatable, save :: MoL_vars(:)
   integer, save :: MoL_nvars = 0

   ! For convenience
   integer, allocatable, save :: MoL_unk_mask(:), MoL_scratch_mask(:)

end module MoL_variables
