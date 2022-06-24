!!****if* source/numericalTools/MoL/MoLMain/ERK/ml_calcRHS
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
!!      ml_calcRHS
!!
!!  SYNOPSIS
!!
!!      call ml_calcRHS(integer, intent(in) :: rhsType
!!                      integer, intent(in) :: rhsStruct
!!                      real,    intent(in) :: t)
!!
!!  DESCRIPTION
!!
!!      Calculate specified RHS type and store in  requested struct
!!
!!  ARGUMENTS
!!
!!      rhsType   : The type of RHS, one of:
!!                  - MOL_RHS_EXPLICIT
!!                  - MOL_RHS_IMPLICIT
!!                  - MOL_RHS_FAST
!!      rhsStruct : MoL memory data-struct to store RHS in
!!      t         : The time of the RHS is to be evaluated at
!!
!!***
subroutine ml_calcRHS(rhsType, rhsStruct, t)
    use MoL_functions

    use ml_memInterface, only: ml_memGetDataPtr, ml_memReleaseDataPtr, ml_memZero

    use Grid_iterator,  only: Grid_iterator_t
    use Grid_tile,      only: Grid_tile_t
    use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    integer, intent(in) :: rhsType, rhsStruct
    real,    intent(in) :: t

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t) :: tileDesc

    real, dimension(:,:,:,:), pointer :: U, rhs

    ! Zero-out RHS memory
    call ml_memZero(rhsStruct)

    call Grid_getTileIterator(itor, LEAF)

    TileLoop: do
        if (.not. itor%isValid()) exit TileLoop

        call itor%currentTile(tileDesc)

        call ml_memGetDataPtr(tileDesc, U, MOL_EVOLVED)
        call ml_memGetDataPtr(tileDesc, rhs, rhsStruct)

        ! All RHS types are calculated and stored-in/added-to one location
        call MoL_rhsE(tileDesc, rhs, U, t)
        call MoL_rhsI(tileDesc, rhs, U, t)
        call MoL_rhsF(tileDesc, rhs, U, t)

        call ml_memReleaseDataPtr(tileDesc, rhs, rhsStruct)
        call ml_memReleaseDataPtr(tileDesc, U, MOL_EVOLVED)

        call itor%next()
    end do TileLoop

    call Grid_releaseTileIterator(itor)
end subroutine ml_calcRHS
