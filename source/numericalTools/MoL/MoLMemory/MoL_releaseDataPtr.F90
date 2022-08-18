!!****if* source/numericalTools/MoL/MoLMemory/MoL_releaseDataPtr
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
!!  NAME
!!
!!      MoL_releaseDataPtr
!!
!!  SYNOPSIS
!!
!!      call MoL_releaseDataPtr(class(Grid_tile_t), intent(in) :: tileDesc
!!                              real, pointer                  :: dataPtr
!!                              integer, intent(in)            :: dataStruct)
!!
!!  DESCRIPTION
!!
!!      Release pointer to the requested data struct for the provided tile
!!
!!      Valid data structs include (defined in MoL.h):
!!          - MOL_EVOLVED : Evolved variables in UNK
!!          - MOL_INITIAL : Copy of the evolved variables at the start of a timestep
!!          - MOL_RHS     : The currently-being-calculated RHS terms
!!          - other       : Each integrator may specify some additional number of
!!                          of scratch-memory for intermediate stages/RHS terms
!!
!!  ARGUMENTS
!!
!!      tileDesc   : Grid tile-descriptor
!!      dataPtr    : Pointer to release
!!      dataStruct : Which data struct
!!
!!***
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
