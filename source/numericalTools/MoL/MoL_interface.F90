!!****h* source/numericalTools/MoL/MoL_interface
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
!!      MoL_interface
!!
!!  SYNOPSIS
!!
!!      use MoL_interface
!!
!!  DESCRIPTION
!!
!!      This is the header file for the method of lines time integration unit
!!      that defines its public interfaces.
!!
!!***

module MoL_interface

    use Grid_tile, only: Grid_tile_t

    implicit none

    !! ================================= !!
    !!  Initialization and finalization  !!
    !! ================================= !!

    interface
        subroutine MoL_init
        end subroutine MoL_init
    end interface

    interface
        subroutine MoL_finalize
        end subroutine MoL_finalize
    end interface

    !! ======================= !!
    !!  Variable registration  !!
    !! ======================= !!

    interface
        subroutine MoL_registerVariable(name, evolIndex, rhsIndex)
            character(len=*), intent(in) :: name
            integer, intent(in) :: evolIndex
            integer, intent(out) :: rhsIndex
        end subroutine MoL_registerVariable
    end interface

    !! ===================== !!
    !!  Advance a time step  !!
    !! ===================== !!

    interface
        subroutine MoL_advance(t, dt)
            real, intent(in) :: t, dt
        end subroutine MoL_advance
    end interface

    !! ================= !!
    !!  MoL data access  !!
    !! ================= !!

    interface
        subroutine MoL_getDataPtr(tileDesc, dataPtr, dataStruct)
            import :: Grid_tile_t
            class(Grid_tile_t), intent(inout) :: tileDesc
            real, dimension(:,:,:,:), pointer :: dataPtr
            integer, intent(in) :: dataStruct
        end subroutine MoL_getDataPtr
    end interface

    interface
        subroutine MoL_releaseDataPtr(tileDesc, dataPtr, dataStruct)
            import :: Grid_tile_t
            class(Grid_tile_t), intent(inout) :: tileDesc
            real, dimension(:,:,:,:), pointer :: dataPtr
            integer, intent(in) :: dataStruct
        end subroutine MoL_releaseDataPtr
    end interface

    !! =================================================== !!
    !!  Call after regrid to update MoL memory structures  !!
    !! =================================================== !!
    interface
        subroutine MoL_regrid
        end subroutine MoL_regrid
    end interface

end module MoL_interface
