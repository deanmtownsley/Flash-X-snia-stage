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
!! @brief MoL_releaseDataPtr stub

!> @ingroup MoL
!!
!! @brief Release a pointer to a MoL data structure for the current tile
!!
!! @details
!! @anchor MoL_releaseDataPtr_stub
!!
!! For compatibility with the reference counting in Grid, all calls to this
!! procedure must be matched by calls to @ref mol_getdataptr
!!
!! @todo This is intended as a temporary measure until a more suitable
!!       solution for MoL's scratch memory is decided
!!
!! @param tileDesc   Descriptor for the current grid tile
!! @param dataPtr    Pointer that targets the specified data structure
!!                   (null on return)
!! @param dataStruct Identifier of the MoL data structure
subroutine MoL_releaseDataPtr(tileDesc, dataPtr, dataStruct)
   use Grid_tile, only: Grid_tile_t

   implicit none

   class(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: dataPtr(:, :, :, :)
   integer, intent(in) :: dataStruct

   nullify (dataPtr)

   return
end subroutine MoL_releaseDataPtr
