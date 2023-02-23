!> @file
!! @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Simulation.h"

!> @ingroup GridMilhoja
!! @stubref{Grid_fillGuardCells}
!!
!! @brief Concrete implementation of Grid_fillGuardCells
!!
!! @attention
!! Only partial functionality has been implemented so far.  This routine
!! aborts if calling code attempts to use unimplemented functionality.
!!
!! @todo Code up full implementation
!! @todo gr_setGcFillNLayers appears to be Paramesh-specific.  Can we get rid
!!       of it?
subroutine Grid_fillGuardCells(gridDataStruct, idir, &
                               minLayers, &
                               eosMode, doEos, &
                               maskSize, mask, makeMaskConsistent, doLogMask, &
                               selectBlockType, &
                               unitReadsMeshDataOnly)
    use milhoja_types_mod,   ONLY : MILHOJA_INT
    use milhoja_grid_mod,    ONLY : milhoja_grid_fillGuardCells

    use Driver_interface,    ONLY : Driver_abort
    use Grid_interface,      ONLY : Grid_getTileIterator, &
                                    Grid_releaseTileIterator
    use Grid_data,           ONLY : gr_justExchangedGC, &
                                    gr_eosMode, &
                                    gr_gcellsUpToDate, &
                                    gr_meshMe
    use Grid_iterator,       ONLY : Grid_iterator_t
    use Grid_tile,           ONLY : Grid_tile_t
    use gr_interface,        ONLY : gr_setGcFillNLayers
    use gr_milhojaInterface, ONLY : gr_checkMilhojaError
    use Eos_interface,       ONLY : Eos_guardCells

    integer, intent(IN)           :: gridDataStruct
    integer, intent(IN)           :: idir
    integer, intent(IN), optional :: minLayers
    integer, intent(IN), optional :: eosMode
    logical, intent(IN), optional :: doEos
    integer, intent(IN), optional :: maskSize
    logical, intent(IN), optional :: mask(:)
    logical, intent(IN), optional :: makeMaskConsistent
    logical, intent(IN), optional :: doLogMask
    integer, intent(IN), optional :: selectBlockType
    logical, intent(IN), optional :: unitReadsMeshDataOnly

    type(Grid_iterator_t)         :: itor
    type(Grid_tile_t)             :: tileDesc
    real,                 pointer :: solnData(:,:,:,:)

    logical, save :: isFirstCall = .TRUE.

    integer :: gcEosMode
    logical :: needEos
    integer :: layers(1:MDIM)
    integer :: returnLayers(1:MDIM)

    integer(MILHOJA_INT) :: MH_err

    NULLIFY(solnData)

    if (gridDataStruct /= CENTER) then
        call Driver_abort("[Grid_fillGuardCells] Unsupported gridDataStruct")
    else if (idir /= ALLDIR) then
        call Driver_abort("[Grid_fillGuardCells] idir must be ALLDIR")
    else if (present(unitReadsMeshDataOnly)) then
        write(*,*) "WARNING: Ignoring unitReadsMeshDataOnly as not yet implemented"
    end if

    if (isFirstCall .AND. (gr_meshMe == MASTER_PE)) then
        if      (present(maskSize)) then
            write(*,*) "WARNING: Ignoring maskSize as not yet implemented"
        else if (present(mask)) then
            write(*,*) "WARNING: Ignoring request to mask as not yet implemented"
        else if (present(makeMaskConsistent)) then
            write(*,*) "WARNING: Ignoring makeMaskConsistent as not yet implemented"
        else if (present(doLogMask)) then
            write(*,*) "WARNING: Ignoring doLogMask as not yet implemented"
        else if (present(selectBlockType)) then
            write(*,*) "WARNING: Ignoring selectBlockType as not yet implemented"
        end if
    end if

    if (present(eosMode)) then
        gcEosMode=eosMode
    else
        gcEosMode=gr_eosMode
    end if

    if (present(doEos)) then
        needEos = doEos
    else
        needEos = .FALSE.
    end if

    !----------------------------------------------------------------
    ! Figure out nlayers arguments to amr_guardcell based on our arguments
    CALL gr_setGcFillNLayers(layers, idir, NGUARD, minLayers, returnLayers)

    !!!!! Cell-centered data first
    if ((gridDataStruct == CENTER) .OR. (gridDataStruct == CENTER_FACES)) then
        CALL milhoja_grid_fillGuardCells(MH_err) 
        CALL gr_checkMilhojaError("Grid_fillGuardCells", MH_err)

        !!!!! FINALIZE CELL-CENTERED DATA
        ! Clean data to account for possible unphysical values caused by
        ! interpolation, revert to primitive form if needed, and
        ! run EoS if needed
        if (needEos) then
            CALL Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
            do while (itor%isValid())
                CALL itor%currentTile(tileDesc)
                CALL tileDesc%getDataPtr(solnData, CENTER)
                
                ! This call disallows the use of tiling
                CALL Eos_guardCells(gcEosMode, solnData, corners=.true., &
                                    layers=returnLayers, &
                                    blockDesc=tileDesc)
                
                CALL tileDesc%releaseDataPtr(solnData, CENTER)
                CALL itor%next()
            end do
            CALL Grid_releaseTileIterator(itor)
        end if
    end if
    
    gr_justExchangedGC = .TRUE.
    gr_gcellsUpToDate = .FALSE.

    isFirstCall = .FALSE.
end subroutine Grid_fillGuardCells

