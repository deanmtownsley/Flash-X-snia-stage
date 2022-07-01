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
!! This is a Milhoja-specific implementation of this routine.  Please refer
!! to the documentation in this routine's stub for general interface information.
!!
!! This routine is needed in order to compile test Milhoja simulations.  However,
!! it aborts if called since it is not yet implemented.
!!
!! @todo Implement.  Should it be in a dedicated file or should changes be made
!!       to the version one folder up?
!! @todo Does this file need the REORDER directive?
subroutine gr_bcPutRegion(gridDataStruct, axis, endPoints, regionSize, mask, &
                          region, tileDesc, idest)
  
    use Driver_interface, ONLY : Driver_abort
    use Grid_tile,        ONLY : Grid_tile_t

    implicit none

    integer,           intent(IN) :: gridDataStruct
    integer,           intent(IN) :: axis
    integer,           intent(IN) :: endPoints(LOW:HIGH, 1:MDIM)
    integer,           intent(IN) :: regionSize(1:REGION_DIM)
    logical,           intent(IN) :: mask(1:regionSize(STRUCTSIZE))
    real,              intent(IN) :: region(1:regionSize(BC_DIR),     &
                                            1:regionSize(SECOND_DIR), &
                                            1:regionSize(THIRD_DIR),  &
                                            1:regionSize(STRUCTSIZE))
    type(Grid_tile_t), intent(IN) :: tileDesc
    integer,           intent(IN) :: idest

    CALL Driver_abort("[gr_bcPutRegion] Not implemented yet")
end subroutine gr_bcPutRegion

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
!! This is a Milhoja-specific implementation of this routine.  Please refer
!! to the documentation in this routine's stub for general interface information.
!!
!! This routine is needed in order to compile test Milhoja simulations.  However,
!! it aborts if called since it is not yet implemented.
!!
!! @todo Implement
subroutine gr_bcPutRegionsMixedGds(gridDataStruct, axis, secondDir, thirdDir, endPoints, &
                                   regionSize, &
                                   regionC, regionFN, regionFT1, regionFT2, &
                                   tileDesc, idest)
  
    use Driver_interface, ONLY : Driver_abort
    use Grid_tile,        ONLY : Grid_tile_t

    implicit none

    integer,           intent(IN)         :: gridDataStruct
    integer,           intent(IN)         :: axis
    integer,           intent(IN)         :: secondDir
    integer,           intent(IN)         :: thirdDir
    integer,           intent(IN)         :: endPoints(LOW:HIGH, 1:MDIM)
    integer,           intent(IN)         :: regionSize(1:REGION_DIM)
    real,                         pointer :: regionC(:,:,:,:)
    real,                         pointer :: regionFN(:,:,:,:)
    real,                         pointer :: regionFT1(:,:,:,:)
    real,                         pointer :: regionFT2(:,:,:,:)
    type(Grid_tile_t), intent(IN)         :: tileDesc
    integer,           intent(IN)         :: idest

    NULLIFY(regionC)
    NULLIFY(regionFN)
    NULLIFY(regionFT1)
    NULLIFY(regionFT2)
    CALL Driver_abort("[gr_bcPutRegionsMixedGds] Not implemented yet")
end subroutine gr_bcPutRegionsMixedGds

