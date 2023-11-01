!!****if* source/physics/IncompNS/IncompNSMain/IncompNS_reInitGridVars
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
!! NAME
!!
!!  IncompNS_reInitGridVars
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!***

!!REORDER(4): face[xyz]Data
!!REORDER(4): solnData

#include "constants.h"
#include "IncompNS.h"
#include "Simulation.h"

subroutine IncompNS_reInitGridVars(tileDesc)

   use Grid_tile, ONLY: Grid_tile_t
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use IncompNS_data

   !------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   integer :: i, j, k
   !------------------------------------------------------------------------------------------
   nullify (solnData, facexData, faceyData, facezData)

   call Timers_start("IncompNS_reInitGridVars")

   ! DEVNOTE (10/24/2023):
   ! Perform allocations based on the unit configuration
   ! cell-centered is only to be initialized when using
   ! IncompNS in variable density mode. Face-centered
   ! in z-direction only required when NDIM == MDIM
   !
   ! Loop over internal cells and set initialization for each data structure
   ! This loop is expensive than it needs to be. Need to optimize or remove
   ! initialization altogether by modifying the algorithm elsewhere in the
   ! code. Currently these initialization are needed for computing the correct
   ! solution, and they will be needed in the future as well but their placement
   ! can be improved.
   !
   ! Start with cell-centered data. Only need for variable density configuration
#ifdef INCOMPNS_VARDENS
   call tileDesc%getDataPtr(solnData, CENTER)
   do k = tileDesc%blkLimitsGC(LOW, KAXIS), tileDesc%blkLimitsGC(HIGH, KAXIS)
      do j = tileDesc%blkLimitsGC(LOW, JAXIS), tileDesc%blkLimitsGC(HIGH, JAXIS)
         do i = tileDesc%blkLimitsGC(LOW, IAXIS), tileDesc%blkLimitsGC(HIGH, IAXIS)
            solnData(PRES_VAR, i, j, k) = 0.
            solnData(RHOC_VAR, i, j, k) = 1.
            solnData(VISC_VAR, i, j, k) = 1.
         end do
      end do
   end do
   call tileDesc%releaseDataPtr(solnData, CENTER)
#endif

   ! face-centered data in x direction, loop over an extra cell in the
   ! corresponding dimension. Also add directive for variable density configuration
   call tileDesc%getDataPtr(facexData, FACEX)
   do k = tileDesc%blkLimitsGC(LOW, KAXIS), tileDesc%blkLimitsGC(HIGH, KAXIS)
      do j = tileDesc%blkLimitsGC(LOW, JAXIS), tileDesc%blkLimitsGC(HIGH, JAXIS)
         do i = tileDesc%blkLimitsGC(LOW, IAXIS), tileDesc%blkLimitsGC(HIGH, IAXIS) + 1
            facexData(HVN0_FACE_VAR, i, j, k) = 0.
#ifdef INCOMPNS_VARDENS
            facexData(RHOF_FACE_VAR, i, j, k) = 1.
            facexData(SIGM_FACE_VAR, i, j, k) = 0.
#endif
            facexData(VFRC_FACE_VAR, i, j, k) = 0.
         end do
      end do
   end do
   call tileDesc%releaseDataPtr(facexData, FACEX)

   ! face-centered data in y direction, loop over extra cell for index j
   call tileDesc%getDataPtr(faceyData, FACEY)
   do k = tileDesc%blkLimitsGC(LOW, KAXIS), tileDesc%blkLimitsGC(HIGH, KAXIS)
      do j = tileDesc%blkLimitsGC(LOW, JAXIS), tileDesc%blkLimitsGC(HIGH, JAXIS) + 1
         do i = tileDesc%blkLimitsGC(LOW, IAXIS), tileDesc%blkLimitsGC(HIGH, IAXIS)
            faceyData(HVN0_FACE_VAR, i, j, k) = 0.
#ifdef INCOMPNS_VARDENS
            faceyData(RHOF_FACE_VAR, i, j, k) = 1.
            faceyData(SIGM_FACE_VAR, i, j, k) = 0.
#endif
            faceyData(VFRC_FACE_VAR, i, j, k) = 0.
         end do
      end do
   end do
   call tileDesc%releaseDataPtr(faceyData, FACEY)

   ! face-centered data in z direction, peform only when NDIM == MDIM
   ! and loop over extra cell for index k
#if NDIM == MDIM
   call tileDesc%getDataPtr(facezData, FACEZ)
   do k = tileDesc%blkLimitsGC(LOW, KAXIS), tileDesc%blkLimitsGC(HIGH, KAXIS) + 1
      do j = tileDesc%blkLimitsGC(LOW, JAXIS), tileDesc%blkLimitsGC(HIGH, JAXIS)
         do i = tileDesc%blkLimitsGC(LOW, IAXIS), tileDesc%blkLimitsGC(HIGH, IAXIS)
            facezData(HVN0_FACE_VAR, i, j, k) = 0.
#ifdef INCOMPNS_VARDENS
            facezData(RHOF_FACE_VAR, i, j, k) = 1.
            facezData(SIGM_FACE_VAR, i, j, k) = 0.
#endif
            facezData(VFRC_FACE_VAR, i, j, k) = 0.
         end do
      end do
   end do
   call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

   call Timers_stop("IncompNS_reInitGridVars")

end subroutine IncompNS_reInitGridVars
