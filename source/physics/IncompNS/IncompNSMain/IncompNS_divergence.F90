!!****if* source/physics/IncompNS/IncompNSMain/IncompNS_divergence
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
!!***
!!REORDER(4): face[xyz]Data
!!REORDER(4): solnData

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

subroutine IncompNS_divergence(tileDesc)

   use Grid_tile, ONLY: Grid_tile_t
   use ins_interface, ONLY: ins_divergence
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use IncompNS_data

!-----------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   integer, dimension(2, MDIM) :: stnLimits = 1
#if NDIM < MDIM
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData
   real, dimension(NFACE_VARS, 1, 1, 1) :: facezData
#else
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
#endif
   real del(MDIM)
   integer NStep
!----------------------------------------------------------------------------------------
#if NDIM < MDIM
   nullify (solnData, facexData, faceyData)
#else
   nullify (solnData, facexData, faceyData, facezData)
#endif

   call Timers_start("IncompNS_divergence")

   stnLimits(LOW, 1:NDIM) = tileDesc%limits(LOW, 1:NDIM) - tileDesc%blkLimitsGC(LOW, 1:NDIM) + 1
   stnLimits(HIGH, 1:NDIM) = tileDesc%limits(HIGH, 1:NDIM) - tileDesc%blkLimitsGC(LOW, 1:NDIM) + 1

   call tileDesc%deltas(del)

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%getDataPtr(facexData, FACEX)
   call tileDesc%getDataPtr(faceyData, FACEY)

#if NDIM == 3
   call tileDesc%getDataPtr(facezData, FACEZ)
#endif

   ! compute divergence of intermediate velocities
   call ins_divergence(facexData(VELC_FACE_VAR, :, :, :), &
                       faceyData(VELC_FACE_VAR, :, :, :), &
                       facezData(VELC_FACE_VAR, :, :, :), &
                       stnLimits(LOW, IAXIS), stnLimits(HIGH, IAXIS), &
                       stnLimits(LOW, JAXIS), stnLimits(HIGH, JAXIS), &
                       stnLimits(LOW, KAXIS), stnLimits(HIGH, KAXIS), &
                       del(DIR_X), del(DIR_Y), del(DIR_Z), &
                       solnData(DUST_VAR, :, :, :))

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)

#if NDIM ==3
   call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

   call Timers_stop("IncompNS_divergence")

   return
end subroutine IncompNS_divergence
