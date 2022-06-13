!!****f* source/numericalTools/MoL/MoL_registerVariable
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
!!      MoL_registerVariable
!!
!!  SYNOPSIS
!!
!!      call MoL_registerVariable(character, intent(in)  :: name(:).
!!                                integer,   intent(in)  :: evolIndex
!!                                integer,   intent(out) :: rhsIndex)
!!
!!  DESCRIPTION 
!!
!!      Register an evolved variable/equation
!!
!! ARGUMENTS
!!
!!      name      : Name of the evolved variable
!!      evolIndex : Index of the evolved variable in UNK
!!      rhsIndex  : Index of the evolved variable in RHS memory
!!
!!***
subroutine MoL_registerVariable(name, evolIndex, rhsIndex)
    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: evolIndex
    integer, intent(out) :: rhsIndex

    rhsIndex = -1

    return
end subroutine MoL_registerVariable
