!!****if* source/physics/ImBound/ImBoundMain/ImBound_reInitGridVars
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
!!
!!
!!
!! ImBound_reInitGridVars
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!***
!!REORDER(4): solnData

#include "constants.h"
#include "ImBound.h"
#include "Simulation.h"

subroutine ImBound_reInitGridVars(tileDesc)

   use Grid_tile, ONLY: Grid_tile_t
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use ImBound_data, ONLY: ib_meshMe

   !------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   real, pointer, dimension(:, :, :, :) :: solnData
   integer :: i, j, k
   !------------------------------------------------------------------------------------------
   nullify (solnData)

   call Timers_start("ImBound_reInitGridVars")

   call tileDesc%getDataPtr(solnData, CENTER)

   do k = tileDesc%blkLimitsGC(LOW, KAXIS), tileDesc%blkLimitsGC(HIGH, KAXIS)
      do j = tileDesc%blkLimitsGC(LOW, JAXIS), tileDesc%blkLimitsGC(HIGH, JAXIS)
         do i = tileDesc%blkLimitsGC(LOW, IAXIS), tileDesc%blkLimitsGC(HIGH, IAXIS)

            ! DEVNOTE (10/24/2023):
            ! Updating intializations for relevant
            ! physics units. See IncompNS_reInitGridVars
            ! for explanation
            solnData(HLN0_VAR, i, j, k) = 0.
         end do
      end do
   end do

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call Timers_stop("ImBound_reInitGridVars")

end subroutine ImBound_reInitGridVars
