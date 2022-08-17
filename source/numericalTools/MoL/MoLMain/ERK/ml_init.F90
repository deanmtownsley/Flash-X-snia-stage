!!****if* source/numericalTools/MoL/MoLMain/ERK/ml_init
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
!!      ml_init
!!
!!  SYNOPSIS
!!
!!      call ml_init()
!!
!!  DESCRIPTION
!!
!!      Initialize a method of lines unit implementation
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine ml_init()
   use MoL_data, only: MoL_nscratch
   use erk_data
   use erk_tableau, only: erk_tableau_init

   use RuntimeParameters_interface, only: RuntimeParameters_get

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

   implicit none

   character(len=MAX_STRING_LENGTH) :: erk_method_str

   integer :: i

   call RuntimeParameters_get("erk_method", erk_method_str)
   erk_method = trim(erk_method_str)

   call erk_tableau_init

   ! Setup RHS indexing
   allocate (erk_K(erk_stages))

   erk_K = (/(MOL_RHS + i, i=0, erk_stages - 1)/)

   MoL_nscratch = erk_stages - 1
end subroutine ml_init
