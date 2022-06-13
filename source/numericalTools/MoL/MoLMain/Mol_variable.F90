!!****ih* source/numericalTools/MoL/MoLMain/MoL_variable
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
!!      MoL_variable
!!
!!  SYNOPSIS
!!
!!      use MoL_variable
!!
!!  DESCRIPTION
!!
!!      Variable-descriptor type used internally by MoL
!!
!!***

module MoL_variable

    implicit none

    type :: MoL_variable_t
        character(len=:), allocatable :: name
        integer :: evolIndex, rhsIndex
    end type MoL_variable_t

end module MoL_variable