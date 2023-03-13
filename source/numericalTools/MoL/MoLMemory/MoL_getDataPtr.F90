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
!! @brief MoL_getDataPtr implementation

!> @ingroup MoLMemory
!!
!! @brief Implements MoL_getDataPtr
!!
!! @details
!! This implementation works only with UG and Paramesh due to the assumed
!! form of the tile-descriptor
!!
!! @todo Implement reference-counting and strictly require that `dataPtr`
!!       is null on entry
!!
!! @stubref{MoL_getDataPtr}

!!REORDER(4): dataPtr
subroutine MoL_getDataPtr(tileDesc, dataPtr, dataStruct)
   use ml_memData, only: ml_scratch_data, ml_activeRHS
   use ml_interface, only: ml_error

   use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

   implicit none

   class(Grid_tile_t), intent(in) :: tileDesc
   real, dimension(:, :, :, :), pointer :: dataPtr
   integer, intent(in) :: dataStruct

   integer :: ind

   if (dataStruct .lt. 0) call ml_error("Unsupported data struct requested")

   if (dataStruct .eq. MOL_EVOLVED) then
      if (associated(dataPtr)) call tileDesc%releaseDataPtr(dataPtr, CENTER)

      ! Grab UNK pointer and bail
      call tileDesc%getDataPtr(dataPtr, CENTER)
   else
      ! Determine if MOL_RHS maps to a specific stage/type of RHS
      if ((dataStruct .eq. MOL_RHS) .and. (ml_activeRHS .ne. MOL_INVALID)) then
         ind = ml_activeRHS
      else
         ind = dataStruct
      end if
      if (associated(dataPtr)) nullify (dataPtr)

      ! Grid_tile_t uses `id` to reference the block
      associate (lo => tileDesc%limits(LOW, :))
         dataPtr(1:, lo(IAXIS):, lo(JAXIS):, lo(KAXIS):) &
            => ml_scratch_data(:, :, :, :, tileDesc%id, ind)
      end associate
   end if
end subroutine MoL_getDataPtr
