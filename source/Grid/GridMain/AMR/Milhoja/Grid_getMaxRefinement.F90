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
!! @todo Confirm that what is coded is correct and then code up the 
!!       remaining functionality.
subroutine Grid_getMaxRefinement(maxRefinement, mode, scope, inputComm)
    use milhoja_types_mod, ONLY : MILHOJA_INT
    use milhoja_grid_mod,  ONLY : milhoja_grid_getCurrentFinestLevel

    use Grid_data,           ONLY : gr_lRefineMax
    use gr_milhojaInterface, ONLY : gr_checkMilhojaError
    use Driver_interface,    ONLY : Driver_abort
    
    implicit none
    
    integer, intent(OUT)          :: maxRefinement
    integer, intent(IN), OPTIONAL :: mode
    integer, intent(IN), OPTIONAL :: scope
    integer, intent(IN), OPTIONAL :: inputComm

    integer :: myMode

    integer(MILHOJA_INT) :: MH_level
    integer(MILHOJA_INT) :: MH_ierr

    maxRefinement = -1

    if (present(mode)) then
        myMode = mode
    else
        myMode = 3
    end if

    if (present(scope)) then
        CALL Driver_abort("[Grid_getMaxRefinement] scope not coded yet")
    end if

    if (present(inputComm)) then
        CALL Driver_abort("[Grid_getMaxRefinement] inputComm not coded yet")
    end if

    if      (myMode == 1) then
        maxRefinement = gr_lRefineMax
    else if (myMode == 4) then
        CALL milhoja_grid_getCurrentFinestLevel(MH_level, MH_ierr)
        CALL gr_checkMilhojaError("Grid_getMaxRefinement", MH_ierr)
        maxRefinement = INT(MH_level)
    else
        CALL Driver_abort("[Grid_getMaxRefinement] mode not implemented")
    end if
end subroutine Grid_getMaxRefinement

