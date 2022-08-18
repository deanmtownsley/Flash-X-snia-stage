!!****if* source/numericalTools/MoL/MoLMemory/MoL_getDataPtr
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
!!      MoL_getDataPtr
!!
!!  SYNOPSIS
!!
!!      call MoL_getDataPtr(class(Grid_tile_t), intent(in) :: tileDesc
!!                          real, pointer                  :: dataPtr
!!                          integer, intent(in)            :: dataStruct)
!!
!!  DESCRIPTION
!!
!!      Obtain pointer to the requested data struct for the provided tile
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
!!      dataPtr    : Pointer to set
!!      dataStruct : Which data struct
!!
!!***
!!REORDER(4): dataPtr
subroutine MoL_getDataPtr(tileDesc, dataPtr, dataStruct)
   use ml_memData, only: scratch_data, ml_activeRHS
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
      if ((dataStruct .eq. MOL_RHS) .and. (ml_activeRHS .ne. MOL_INVALID)) then
         ind = ml_activeRHS
      else
         ind = dataStruct
      end if
      if (associated(dataPtr)) nullify (dataPtr)

      ! Grid_tile_t uses `id` to reference the block
      associate (lo => tileDesc%limits(LOW, :))
         dataPtr(1:, lo(IAXIS):, lo(JAXIS):, lo(KAXIS):) &
            => scratch_data(:, :, :, :, tileDesc%id, ind)
      end associate
   end if
end subroutine MoL_getDataPtr
