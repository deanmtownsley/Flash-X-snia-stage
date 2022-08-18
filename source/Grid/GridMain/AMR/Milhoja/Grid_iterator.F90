#include "constants.h"
#include "Simulation.h"

!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! This is a Milhoja-specific implementation of the Grid_iterator_t class. While
!! its public interface must match that of the Grid_iterator_t classes defined
!! by other Grid implementations, its implementation can be significantly
!! different.  Note that there is no mechanism to enforce this interface
!! requirement.
!!
!! This implementation is presently limited to pseudo-UG only operation.
!!
!! @todo Code up full implementation
module Grid_iterator
    use iso_c_binding,     ONLY : C_PTR, &
                                  C_NULL_PTR

    use milhoja_types_mod, ONLY : MILHOJA_INT

    implicit none

    private

    public :: build_iterator
    public :: destroy_iterator

    !> Instances of this class can be used to iterate over and access a subset
    !! of tiles in the grid.
    type, public :: Grid_iterator_t
        type(C_PTR), private :: MH_itor_Cptr = C_NULL_PTR
    contains
        procedure, public :: isValid
        procedure, public :: next
        procedure, public :: currentTile
    end type Grid_iterator_t

    !!!!!----- INTERFACES TO C-LINKAGE C++ FUNCTIONS
    ! The C-to-Fortran interoperability layer
    interface
        !> Fortran interface on routine in C interface of same name.
        function milhoja_itor_build_C(C_itor_ptr) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            implicit none
            type(C_PTR),         intent(OUT) :: C_itor_ptr
            integer(MILHOJA_INT)             :: C_ierr
        end function milhoja_itor_build_C

        !> Fortran interface on routine in C interface of same name.
        function milhoja_itor_destroy_C(C_itor_ptr) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            implicit none
            type(C_PTR),         intent(IN), value :: C_itor_ptr
            integer(MILHOJA_INT)                   :: C_ierr
        end function milhoja_itor_destroy_C

        !> Fortran interface on routine in C interface of same name.
        function milhoja_itor_is_valid_C(C_itor_ptr, C_isValid) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            implicit none
            type(C_PTR),         intent(IN), value :: C_itor_ptr
            type(C_PTR),         intent(IN), value :: C_isValid
            integer(MILHOJA_INT)                   :: C_ierr
        end function milhoja_itor_is_valid_C

        !> Fortran interface on routine in C interface of same name.
        function milhoja_itor_next_C(C_itor_ptr) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            implicit none
            type(C_PTR),         intent(IN), value :: C_itor_ptr
            integer(MILHOJA_INT)                   :: C_ierr
        end function milhoja_itor_next_C

        !> Fortran interface on routine in C interface of same name.
        function milhoja_itor_acquire_current_tile_C(C_itor_ptr, &
                                                     C_tile_ptr) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            implicit none
            type(C_PTR),         intent(IN), value :: C_itor_ptr
            type(C_PTR),         intent(OUT)       :: C_tile_ptr
            integer(MILHOJA_INT)                   :: C_ierr
        end function milhoja_itor_acquire_current_tile_C

        !> Fortran interface on routine in C interface of same name.
        function milhoja_itor_release_current_tile_C(C_tile_ptr) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            implicit none
            type(C_PTR),         intent(IN), value :: C_tile_ptr
            integer(MILHOJA_INT)                   :: C_ierr
        end function milhoja_itor_release_current_tile_C
    end interface

contains

    !> Build and obtain a tile iterator.  As part of the build, resources are
    !! obtained and the iterator owns these resources.  Therefore, code that
    !! builds an iterator must destroy the iterator using destroy_iterator()
    !! once the iterator is no longer needed.
    !!
    !! It is the responsibility of calling code to ensure that the iterator is
    !! used in common and reasonable ways.  For instance, the iterator should
    !! not be used if the Grid data structures might have been altered by
    !! actions such as regridding after the iterator was acquired.
    !!
    !! @todo Implement full functionality.
    !! @todo This documentation is mostly generic and should be put into a stub
    !!       file.
    !!
    !! @param itor       The iterator
    !! @param nodetype   The types of tiles that should be accessed by the 
    !!                   iterator.  Presently, only LEAF and ALL_BLKS are
    !!                   acceptable, but both access all blocks.
    subroutine build_iterator(itor, nodetype)
        use gr_milhojaInterface, ONLY : gr_checkMilhojaError
        use Driver_interface,    ONLY : Driver_abort

        type(Grid_iterator_t), intent(OUT) :: itor
        integer,               intent(IN)  :: nodetype

        integer(MILHOJA_INT) :: MH_ierr

        if ((nodetype /= LEAF) .AND. (nodetype /= ALL_BLKS)) then
            CALL Driver_abort("[build_iterator] nodetype not supported yet")
        end if

        MH_ierr = milhoja_itor_build_C(itor%MH_itor_Cptr)
        CALL gr_checkMilhojaError("build_iterator", MH_ierr)
    end subroutine build_iterator

    !> Release the resources associated with the given iterator and destroy the
    !! iterator.
    !!
    !! @todo Prepending the subroutine line with IMPURE_ELEMENTAL confuses
    !!       Doxygen.  IMPURE_ELEMENTAL was copied directly from AMReX.  What
    !!       does it do and is it necessary?
    !!
    !! @param itor       The iterator to destroy
!    IMPURE_ELEMENTAL subroutine destroy_iterator(itor)
    subroutine destroy_iterator(itor)
        use gr_milhojaInterface, ONLY : gr_checkMilhojaError

        type(Grid_iterator_t), intent(INOUT) :: itor

        integer(MILHOJA_INT) :: MH_ierr

        MH_ierr = milhoja_itor_destroy_C(itor%MH_itor_Cptr)
        CALL gr_checkMilhojaError("destroy_iterator", MH_ierr)

        itor%MH_itor_Cptr = C_NULL_PTR
    end subroutine destroy_iterator

    !> Determine if the iterator is valid and, therefore, if calling code can
    !! safely call next().
    !!
    !! @returns True if the iterator is valid and can be advanced with next();
    !! False, if the iterator is set to the last tile and should *not* be
    !! advanced with next().
    function isValid(this)
        use iso_c_binding, ONLY : C_LOC

        use gr_milhojaInterface, ONLY : gr_checkMilhojaError

        class(Grid_iterator_t), intent(IN)        :: this
        logical,                           target :: isValid

        type(C_PTR)          :: isValid_Cptr
        integer(MILHOJA_INT) :: MH_ierr

        isValid = .FALSE.
        isValid_Cptr = C_LOC(isValid)

        MH_ierr = milhoja_itor_is_valid_C(this%MH_itor_Cptr, isValid_Cptr)
        CALL gr_checkMilhojaError("Grid_iterator%isValid", MH_ierr)
    end function isValid

    !> Advance the iterator to the next tile.  Refer to the documentation for
    !! isValid() for more information on the proper usage of this routine.
    subroutine next(this)
        use gr_milhojaInterface, ONLY : gr_checkMilhojaError

        class(Grid_iterator_t), intent(INOUT) :: this

        integer(MILHOJA_INT) :: MH_ierr

        MH_ierr = milhoja_itor_next_C(this%MH_itor_Cptr)
        CALL gr_checkMilhojaError("Grid_iterator%next", MH_ierr)
    end subroutine next

    !> Obtain a tile object that can be used by calling code to access the
    !! data and metadata of the tile that the iterator currently indexes.
    !!
    !! @param tileDesc   The tile object
    subroutine currentTile(this, tileDesc)
        use iso_c_binding,       ONLY : C_NULL_PTR
        use Grid_tile,           ONLY : Grid_tile_t, &
                                        Grid_tile_fromMilhojaTilePtr
        use gr_milhojaInterface, ONLY : gr_checkMilhojaError

        class(Grid_iterator_t), intent(IN)  :: this
        type(Grid_tile_t),      intent(OUT) :: tileDesc

        type(C_PTR)          :: MH_tile_Cptr
        integer(MILHOJA_INT) :: MH_ierr

        ! MH_tile_Cptr assumes ownership of the tile and its resources.
        MH_tile_Cptr = C_NULL_PTR
        MH_ierr = milhoja_itor_acquire_current_tile_C(this%MH_itor_Cptr, MH_tile_Cptr)
        CALL gr_checkMilhojaError("Grid_iterator%currentTile", MH_ierr)

        CALL Grid_tile_fromMilhojaTilePtr(MH_tile_Cptr, tileDesc)

        ! Release ownership of the tile and its resources since metadata cached
        ! in tileDesc.
        MH_ierr = milhoja_itor_release_current_tile_C(MH_tile_Cptr)
        CALL gr_checkMilhojaError("Grid_iterator%currentTile", MH_ierr)
    end subroutine currentTile 
 
end module Grid_iterator

