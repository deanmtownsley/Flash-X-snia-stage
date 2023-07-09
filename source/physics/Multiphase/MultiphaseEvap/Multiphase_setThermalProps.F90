!!***if* source/physics/Multiphase/MultiphaseEvap/Multiphase_setThermalProps
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

subroutine Multiphase_setThermalProps(tileDesc)

   use Multiphase_data
   use Stencils_interface, ONLY: Stencils_lsCenterProps
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t

!------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:, :, :, :) :: solnData
   integer :: ierr, i, j, k
   real del(MDIM)
   real minCellDiag

!------------------------------------------------------------------------------------------
   nullify (solnData)

   call Timers_start("Multiphase_setThermalProps")

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%deltas(del)

   solnData(mph_iAlphaCVar, :, :, :) = solnData(SMHV_VAR, :, :, :)*(mph_thcoGas/(mph_rhoGas*mph_CpGas)) + &
                                       (1 - solnData(SMHV_VAR, :, :, :))*solnData(mph_iAlphaCVar, :, :, :)

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)

   call Timers_stop("Multiphase_setThermalProps")

   return

end subroutine Multiphase_setThermalProps
