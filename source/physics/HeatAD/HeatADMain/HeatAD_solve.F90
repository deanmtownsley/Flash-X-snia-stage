!!***if* source/physics/HeatAD/HeatADMain/HeatAD_solve
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
#include "HeatAD.h"
#include "Simulation.h"

subroutine HeatAD_solve(tileDesc, dt)

   use HeatAD_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t
   use Stencils_interface, ONLY: Stencils_integrateEuler, Stencils_integrateAB2

   implicit none
   include"Flashx_mpi.h"
   real, INTENT(IN) :: dt
   type(Grid_tile_t), INTENT(IN) :: tileDesc

!--------------------------------------------------------------------------------------------
   real ::  del(MDIM)
   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:, :, :, :) :: solnData
   real :: diffusion_coeff

!---------------------------------------------------------------------------------------------
   call Timers_start("HeatAD_solve")

   nullify (solnData)

   diffusion_coeff = ht_invReynolds/ht_Prandtl

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%deltas(del)

   call Stencils_integrateAB2(solnData(TEMP_VAR, :, :, :), &
                              solnData(HTN0_VAR, :, :, :), &
                              solnData(HTN1_VAR, :, :, :), &
                              dt, &
                              GRID_ILO, GRID_IHI, &
                              GRID_JLO, GRID_JHI, &
                              GRID_KLO, GRID_KHI, &
                              iSource=solnData(TFRC_VAR, :, :, :))

   solnData(HTN1_VAR, :, :, :) = solnData(HTN0_VAR, :, :, :)

   call tileDesc%releaseDataPtr(solnData, CENTER)

   call Timers_stop("HeatAD_solve")

end subroutine HeatAD_solve
