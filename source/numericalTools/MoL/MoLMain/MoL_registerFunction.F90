!!****if* source/numericalTools/MoL/MoLMain/MoL_registerFunction
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
!!      MoL_registerFunction
!!
!!  SYNOPSIS
!!
!!      call MoL_registerFunction(integer,   intent(in) :: funcType
!!                                procedure,            :: func)
!!
!!  DESCRIPTION 
!!
!!      Register a function with MoL
!!
!! ARGUMENTS
!!
!!      funcType : One of the following function types defined in MoL.h
!!                  MOL_RHS_EXPLICIT      -  RHS for (slow) explicit terms
!!                  MOL_RHS_IMPLICIT      -  RHS for (slow) implicit terms
!!                  MOL_RHS_FAST          -  RHS for (fast) explicit terms
!!                  MOL_IMPLICIT_UPDATE   -  Implicit updates
!!                  MOL_POST_UPDATE       -  Post-update (slow) per-stage
!!                  MOL_POST_UPDATE_FAST  -  Post-update (fast) per-stage
!!      func     : Procedure to register
!!
!!***
subroutine MoL_registerRHS(rhsType, rhsFunc)
    use Grid_tile, only: Grid_tile_t

    use MoL_functions

    use ml_interface, only: ml_error

#include "MoL.h"

    implicit none

    integer, intent(in) :: rhsType

    interface
        subroutine rhsFunc(tileDesc, dy, y, t)
            import :: Grid_tile_t
            class(Grid_tile_t), intent(in) :: tileDesc
            real, dimension(:,:,:,:), pointer :: dy, y
            real, intent(in) :: t
        end subroutine rhsFunc
    end interface

    select case(rhsType)
    case(MOL_RHS_EXPLICIT)
        if (associated(MoL_rhsE)) nullify(MoL_rhsE)
        MoL_rhsE => rhsFunc

    case(MOL_RHS_IMPLICIT)
        if (associated(MoL_rhsI)) nullify(MoL_rhsI)
        MoL_rhsI => rhsFunc

    case(MOL_RHS_FAST)
        if (associated(MoL_rhsF)) nullify(MoL_rhsF)
        MoL_rhsF => rhsFunc

    case default
        call ml_error("Attempting to register unknown RHS function type")
    end select ! rhsType
end subroutine MoL_registerRHS

subroutine MoL_registerUpdate(updateType, updateFunc)
    use MoL_functions

    use ml_interface, only: ml_error

#include "MoL.h"

    implicit none

    integer, intent(in) :: updateType

    interface
        subroutine updateFunc(t, dt)
            real, intent(in) :: t, dt
        end subroutine updateFunc
    end interface

    select case(updateType)
    case(MOL_IMPLICIT_UPDATE)
        if (associated(MoL_implicitUpdate)) nullify(MoL_implicitUpdate)
        MoL_implicitUpdate => updateFunc

    case default
        call ml_error("Attempting to register unknown update function type")
    end select ! updateType
end subroutine MoL_registerUpdate

subroutine MoL_registerPostUpdate(postUpdateType, postUpdateFunc)
    use MoL_functions

    use ml_interface, only: ml_error

#include "MoL.h"

    implicit none

    integer, intent(in) :: postUpdateType

    interface
        subroutine postUpdateFunc(t)
            real, intent(in) :: t
        end subroutine postUpdateFunc
    end interface

    select case(postUpdateType)
    case(MOL_POST_UPDATE)
        if (associated(MoL_postUpdate)) nullify(MoL_postUpdate)
        MoL_postUpdate => postUpdateFunc

    case(MOL_POST_UPDATE_FAST)
        if (associated(MoL_postUpdateFast)) nullify(MoL_postUpdateFast)
        MoL_postUpdateFast => postUpdateFunc

    case default
        call ml_error("Attempting to register unknown post-update function type")
    end select ! postUpdateType
end subroutine MoL_registerPostUpdate
