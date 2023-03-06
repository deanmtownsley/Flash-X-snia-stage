!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief MoL_releaseDataPtr implementation

!> @ingroup MoLMemory
!!
!! @brief Implements MoL_releaseDataPtr
!!
!! @details
!! This implementation will work will all Grid backends
!!
!! @todo Implement reference-counting instead of only nullifying
!!       `dataPtr` for MoL memory targets
!!
!! @stubref{MoL_releaseDataPtr}
subroutine MoL_releaseDataPtr(tileDesc, dataPtr, dataStruct)
   use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

   implicit none

   class(Grid_tile_t), intent(in) :: tileDesc
   real, dimension(:, :, :, :), pointer :: dataPtr
   integer, intent(in) :: dataStruct

   if (dataStruct .eq. MOL_EVOLVED) then
      call tileDesc%releaseDataPtr(dataPtr, CENTER)
   else
      nullify (dataPtr)
   end if
end subroutine MoL_releaseDataPtr
