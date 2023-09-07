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
!! @brief MoL_getDataPtr stub

!> @ingroup MoL
!!
!! @brief Obtain a pointer to a MoL data structure for the current tile
!!
!! @details
!! @anchor MoL_getDataPtr_stub
!!
!! This procedure will associate the provided pointer to the current tile
!! in the requested data structure.  For compatibility with the reference
!! counting in Grid, all calls to this procedure must be matched by calls
!! to @ref mol_releasedataptr
!!
!! Valid data structures include (defined in MoL.h):
!!    - `MOL_EVOLVED` : Evolved variables in UNK
!!    - `MOL_INITIAL` : Copy of the evolved variables at the start of a timestep
!!    - `MOL_RHS`     : The currently-being-calculated RHS terms
!!    - other         : Each integrator may specify some additional number of
!!                      of scratch-memory for intermediate stages/RHS terms
!!
!! @note Requests for `MOL_RHS` will potentially return a pointer to different
!!       data structures during different integration stanges and/or types
!!       of RHS calculations.  This will occur when a user requests the RHS
!!       pointer during a call to one of the RHS procedures if an integration
!!       scheme sets the active RHS data structure for the current stage.  In
!!       these cases a user's code will receive a pointer to the correct
!!       intermediate-stage RHS storage to add its contributions to.  If an
!!       integration scheme does not utilize this feature, the default `MOL_RHS`
!!       memory structure will always be pointed to, and it is assumed that
!!       if an integration scheme requires intermediate-stage storage, a copy
!!       from `MOL_RHS` to the proper storage will be made
!!
!! @note Requests for `MOL_EVOLVED` are forwarded to the provided
!!       tile descriptor
!!
!! @todo This is intended as a temporary measure until a more suitable
!!       solution for MoL's scratch memory is decided
!!
!! @param tileDesc   Descriptor for the current grid tile
!! @param dataPtr    Pointer that will target the current tile in the
!!                   requested data structure
!! @param dataStruct Identifier of the MoL data structure
subroutine MoL_getDataPtr(tileDesc, dataPtr, dataStruct)
   use Grid_tile, only: Grid_tile_t

   implicit none

   class(Grid_tile_t), intent(in) :: tileDesc
   real, pointer :: dataPtr(:, :, :, :)
   integer, intent(in) :: dataStruct

   nullify (dataPtr)

   return
end subroutine MoL_getDataPtr
