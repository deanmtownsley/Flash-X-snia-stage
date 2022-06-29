#include "constants.h"

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
!! A Milhoja-specific implementation of this routine.  Please refer
!! to the documentation in this routine's stub for general interface
!! information.
!!
!! Only partial functionality has been implemented so far.  This routine
!! aborts if calling code attempts to use non-implemented functionality.
!!
!! @todo Code up full implementation
subroutine Grid_getTileIterator(itor, nodetype, level, tiling, tileSize, nthreads)
    use Grid_data,        ONLY : gr_meshMe
    use Grid_iterator,    ONLY : Grid_iterator_t, &
                                 build_iterator
    use Driver_interface, ONLY : Driver_abort
    
    implicit none
    
    type(Grid_iterator_t), intent(OUT)          :: itor
    integer,               intent(IN)           :: nodetype
    integer,               intent(IN), optional :: level
    logical,               intent(IN), optional :: tiling
    integer,               intent(IN), optional :: tileSize(1:MDIM)
    integer,               intent(IN), optional :: nthreads

    logical, save :: isFirstCall = .TRUE.

    if      (present(level)) then
        if (isFirstCall .AND. (gr_meshMe == MASTER_PE)) then
            write(*,*) "WARNING: This is only working for pseudo-UG mode"
        end if
        if (level /= 1) then
            CALL Driver_abort("[Grid_getTileIterator] Invalid level for pseudo-UG")
        end if
    else if (present(nthreads)) then
        CALL Driver_abort("[Grid_getTileIterator] Multithreading not implemented yet")
    end if

    if (present(tiling)) then
        if (tiling .AND. isFirstCall .AND. (gr_meshMe == MASTER_PE)) then
            write(*,*) "WARNING: Ignoring request to tile as this is not yet implemented"
        end if
    end if

    CALL build_iterator(itor, nodetype)

    isFirstCall = .FALSE.
end subroutine Grid_getTileIterator

