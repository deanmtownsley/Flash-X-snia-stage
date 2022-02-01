!!****if* source/physics/Hydro/HydroMain/hy_memReleaseBlkPtr
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
!! NAME
!!  hy_memReleaseBlkPtr
!!
!! SYNOPSIS
!!
!!  call hy_memReleaseBlkPtr(integer(IN)   :: blockId,
!!                     real(pointer) :: dataPtr(:,:,:,:),
!!                     integer(IN),optional :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  Releases a pointer to a block.
!!  
!! ARGUMENTS 
!!
!!  blockId - ID of the block, should be the same ID was used in the
!!            corresponding Grid_getBlkPtr call.
!!  dataPtr - Pointer to be released.
!!  gridDataStruct - an optional argument that designates the type of grid data
!!                   structure to handle (i.e. facevar, unknown, scratch...)
!!
!!
!! NOTES
!!
!!  This implementation actually just nullifies the pointer.
!!
!!***

subroutine hy_memReleaseBlkPtr_desc(tileDesc, dataPtr, gridDataStruct)
  use Grid_tile, ONLY : Grid_tile_t

#include "FortranLangFeatures.fh"

  implicit none

  type(Grid_tile_t), intent(IN) :: tileDesc
  real, POINTER_INTENT_OUT :: dataPtr(:,:,:,:)
  integer,optional, intent(in) :: gridDataStruct

  nullify(dataPtr)

end subroutine hy_memReleaseBlkPtr_desc

