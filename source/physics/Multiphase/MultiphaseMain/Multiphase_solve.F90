!!***if* source/physics/Multiphase/MultiphaseMain/Multiphase_solve
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
!!***
!!REORDER(4): solnData

#include "constants.h"
#include "Multiphase.h"
#include "Simulation.h"

subroutine Multiphase_solve(tileDesc, dt)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t
   use Stencils_interface, ONLY: Stencils_integrateEuler, Stencils_integrateAB2

   implicit none
   include "Flashx_mpi.h"
   !----------Arugments List---------------
   real, INTENT(IN) :: dt
   type(Grid_tile_t), INTENT(IN) :: tileDesc

!-----------------------------------------------------------------------------------------
   integer, dimension(2, MDIM) :: stnLimits = 1
   real, pointer, dimension(:, :, :, :) :: solnData
   real del(MDIM)
!-----------------------------------------------------------------------------------------
   nullify (solnData)

   call Timers_start("Multiphase_solve")

   call tileDesc%getDataPtr(solnData, CENTER)

   stnLimits(LOW, 1:NDIM) = tileDesc%limits(LOW, 1:NDIM) - tileDesc%blkLimitsGC(LOW, 1:NDIM) + 1
   stnLimits(HIGH, 1:NDIM) = tileDesc%limits(HIGH, 1:NDIM) - tileDesc%blkLimitsGC(LOW, 1:NDIM) + 1

   call Stencils_integrateEuler(solnData(DFUN_VAR, :, :, :), &
                                solnData(HDN0_VAR, :, :, :), &
                                dt, &
                                stnLimits(LOW, IAXIS), stnLimits(HIGH, IAXIS), &
                                stnLimits(LOW, JAXIS), stnLimits(HIGH, JAXIS), &
                                stnLimits(LOW, KAXIS), stnLimits(HIGH, KAXIS), &
                                iSource=solnData(DFRC_VAR, :, :, :))

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call Timers_stop("Multiphase_solve")

   return
end subroutine Multiphase_solve
