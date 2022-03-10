!!****if* source/physics/IncompNS/IncompNSMain/constdens/IncompNS_setupPoisson
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

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

subroutine IncompNS_setupPoisson(tileDesc, dt)

   use Grid_tile, ONLY: Grid_tile_t
   use ins_interface, ONLY: ins_setupPoissonRhs_constdens
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use IncompNS_data

   implicit none
   include "Flashx_mpi.h"
   !---Argument List-------
   real, INTENT(IN) :: dt
   type(Grid_tile_t), INTENT(IN) :: tileDesc

!------------------------------------------------------------------------------------------
   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:, :, :, :) :: solnData
   real del(MDIM)
   integer :: NStep
!------------------------------------------------------------------------------------------
   nullify (solnData)

   call Timers_start("IncompNS_setupPoisson")

   !---POISSON RHS:-------------------------------------------------------------------------------------
   blkLimits = tileDesc%limits
   blkLimitsGC = tileDesc%blkLimitsGC

   call tileDesc%deltas(del)
   call tileDesc%getDataPtr(solnData, CENTER)

   ! Poisson RHS source vector
   call ins_setupPoissonRhs_constdens(solnData(DUST_VAR, :, :, :), dt)

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)

   call Timers_stop("IncompNS_setupPoisson")

   return
end subroutine IncompNS_setupPoisson
