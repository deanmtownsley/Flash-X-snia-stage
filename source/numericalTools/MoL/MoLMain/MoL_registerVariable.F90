!!****if* source/numericalTools/MoL/MoLMain/MoL_registerVariable
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
    use MoL_variables
    use ml_interface, only: ml_warn

    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: evolIndex
    integer, intent(out) :: rhsIndex

    type(MoL_variable_t) :: var

    integer :: ivar

    ! If this variable is registered already, output a warning (if verbose) and exit
    varLoop: do ivar = 1, MoL_nvars
        if (MoL_vars(ivar)%evolIndex .eq. evolIndex) then
            call ml_warn("Duplicate variable registration")
            return
        end if
    end do varLoop
    

    MoL_nvars = MoL_nvars + 1

    var%name = trim(name)
    var%evolIndex = evolIndex
    var%rhsIndex = MoL_nvars

    rhsIndex = MoL_nvars

    if (.not. allocated(MoL_vars)) then
        allocate(MoL_vars(1))
        MoL_vars = var
    else
        MoL_vars = (/MoL_vars, var/)
    endif
end subroutine MoL_registerVariable
