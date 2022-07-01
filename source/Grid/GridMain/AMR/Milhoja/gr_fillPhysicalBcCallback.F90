#include "Milhoja.h"
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
!! Given a region of the cell-centered space, this routine will fill those
!! guardcells outside the domain with data that satisfies the boundary conditions
!! of the problem.  It is assumed that the data for the given region already
!! contains data that is valid in the domain.  It is assumed that if a LOW/HIGH
!! pair of faces have periodic BC, then the layers of data outside these faces have
!! also been filled with correct data.
!!
!! This routine does not assume that the given region is a block.  Rather, it
!! could be a region within a block or spanning more than one block.
!!
!! To be meaningful, it is assumed that all regions that span a domain boundary
!! have at least NGUARD layers of interior data along the axis of the
!! boundary.  It is also assumed that these regions do not require filling more
!! than NGUARD layers of guardcells outside the domain boundary.
!!
!! This routine executes the BC fill using the GridBoundaryConditions subunit.
!! In particular, calling code is first given the opportunity to execute the
!! fill via the routine Grid_bcApplyToRegionSpecialized.  If this routine does
!! not handle the fill, then the fill is done via Grid_bcApplyToRegion.
!!
!! Since this is registered with Milhoja as a callback to be used (potentially)
!! with C++ grid backends, it is part of the Fortran/C++ interoperability layer. As
!! a result it is acceptable that it have a C-linkage-compatible interface.
!!
!! NOTE: This routine has been written only as a proof-of-concept of the
!! Fortran/C++ interoperability layer.  It does
!! *not* write any BC data and, therefore, limits the use of Milhoja to only
!! all-periodic BC.
!!
!! NOTE: The interface of this could be changed substantially once this
!! is under active development.
!!
!! @todo We are presently restricted to using cell-centered only.
!! @todo Should lo/hi really be NDIM?  What if we want to use a backend
!!       other than AMReX?
!! @todo What is the variable indexing scheme in Flash-X and Milhoja?
!! @todo How should Milhoja be configured in terms of BCs above NDIM?
!!
!! @param lo_Cptr      Need to figure out what this should be
!! @param hi_Cptr      Need to figure out what this should be
!! @param MH_level     The refinement level of the given region
!! @param MH_startVar  Need to figure out what this should be
!! @param MH_nVars     Need to figure out what this should be
subroutine gr_fillPhysicalBcCallback(lo_Cptr, hi_Cptr, MH_level, &
                                     MH_startVar, MH_nVars) bind(c)
    use iso_c_binding,     ONLY : C_PTR, &
                                  C_F_POINTER

    use milhoja_types_mod, ONLY : MILHOJA_INT

    use Grid_data,         ONLY : gr_domainBC
    use Driver_interface,  ONLY : Driver_abort

    implicit none

    type(C_PTR),          intent(IN), value :: lo_Cptr
    type(C_PTR),          intent(IN), value :: hi_Cptr
    integer(MILHOJA_INT), intent(IN), value :: MH_level
    integer(MILHOJA_INT), intent(IN), value :: MH_startVar
    integer(MILHOJA_INT), intent(IN), value :: MH_nVars

    integer(MILHOJA_INT), pointer :: MH_lo(:)
    integer(MILHOJA_INT), pointer :: MH_hi(:)
    integer                       :: dataShape(1)

    integer :: lo(1:MDIM)
    integer :: hi(1:MDIM)
    integer :: level
    integer :: startVar
    integer :: nVars

    integer :: axis
    integer :: face

    NULLIFY(MH_lo)
    NULLIFY(MH_hi)

    ! Since we should only be running all-periodic BC problems, Milhoja's grid
    ! backend should not even be calling this.
    CALL Driver_abort("[gr_fillPhysicalBcCallback] Periodic-only BC please")

    !!!!!----- CONVERT ALL GIVEN DATA FOR USE BY FLASH-X
    dataShape(1) = MILHOJA_NDIM
    CALL C_F_POINTER(lo_Cptr, MH_lo, shape=dataShape)
    if (.NOT. ASSOCIATED(MH_lo)) then
        CALL Driver_abort("[gr_fillPhysicalBcCallback] MH_lo is NULL")
    end if
    CALL C_F_POINTER(hi_Cptr, MH_hi, shape=dataShape)
    if (.NOT. ASSOCIATED(MH_hi)) then
        CALL Driver_abort("[gr_fillPhysicalBcCallback] MH_hi is NULL")
    end if

    ! Assume that Milhoja is 0-based spatial indices.  Flash-X is 1-based.
    lo(:) = 1
    hi(:) = 1
    do axis = 1, MILHOJA_NDIM
        lo(axis) = INT(MH_lo(axis)) + 1
        hi(axis) = INT(MH_hi(axis)) + 1
    end do
    level    = INT(MH_level) + 1
    startVar = INT(MH_startVar)
    nVars    = INT(MH_nVars)

    write(*,*) "[gr_fillPhysicalBcCallback] ------------------------"
    write(*,*) "[gr_fillPhysicalBcCallback] Region lo   = ", lo
    write(*,*) "[gr_fillPhysicalBcCallback] Region hi   = ", hi
    write(*,*) "[gr_fillPhysicalBcCallback] Level       = ", level
    write(*,*) "[gr_fillPhysicalBcCallback] Start var   = ", startVar
    write(*,*) "[gr_fillPhysicalBcCallback] N Variables = ", nVars
    do     axis = 1, MILHOJA_NDIM
        do face = LOW, HIGH
            write(*,*) "[gr_fillPhysicalBcCallback] Face = ", face, " / Axis = ", axis
            if (gr_domainBC(face, axis) == PERIODIC) then
                write(*,*) "[gr_fillPhysicalBcCallback] Periodic BC"
            else
                write(*,*) "[gr_fillPhysicalBcCallback] Non-periodic BC"
            end if
        end do
    end do
end subroutine gr_fillPhysicalBcCallback

