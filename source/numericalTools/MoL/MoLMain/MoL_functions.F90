!!****if* source/numericalTools/MoL/MoLMain/MoL_functions
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
!!      MoL_functions
!!
!!  SYNOPSIS
!!
!!      use MoL_functions
!!
!!  DESCRIPTION
!!
!!      Function pointers used internally by MoL
!!
!!***

module MoL_functions

    use Grid_tile, only: Grid_tile_t

    implicit none

    abstract interface
        subroutine MoL_rhs_t(tileDesc, dy, y, t)
            import :: Grid_tile_t
            class(Grid_tile_t), intent(in) :: tileDesc
            real, dimension(:,:,:,:), pointer :: dy, y
            real, intent(in) :: t
        end subroutine MoL_rhs_t

        subroutine MoL_implicitUpdate_t(t, dt)
            real, intent(in) :: t, dt
        end subroutine MoL_implicitUpdate_t

        subroutine MoL_postUpdate_t(t)
            real, intent(in) :: t
        end subroutine MoL_postUpdate_t
    end interface

    procedure(MoL_rhs_t), pointer :: MoL_rhsE => null()
    procedure(MoL_rhs_t), pointer :: MoL_rhsI => null()
    procedure(MoL_rhs_t), pointer :: MoL_rhsF => null()

    procedure(MoL_implicitUpdate_t), pointer :: MoL_implicitUpdate => null()

    procedure(MoL_postUpdate_t), pointer :: MoL_postUpdate => null()
    procedure(MoL_postUpdate_t), pointer :: MoL_postUpdateFast => null()

end module MoL_functions
