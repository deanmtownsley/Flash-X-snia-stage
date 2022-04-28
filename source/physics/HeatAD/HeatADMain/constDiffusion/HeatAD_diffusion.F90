!!***if* source/physics/HeatAD/HeatADMain/constDiffusion/HeatAD_diffusion
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

subroutine HeatAD_diffusion(tileDesc)

   use HeatAD_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t
   use Stencils_interface, ONLY: Stencils_diffusion2d, Stencils_diffusion3d

!--------------------------------------------------------------------------------------------
   implicit none
   include"Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   real ::  del(MDIM)
   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:, :, :, :) :: solnData
   real :: diffusion_coeff

!---------------------------------------------------------------------------------------------
   nullify (solnData)

   call Timers_start("HeatAD_diffusion")

   diffusion_coeff = ht_invReynolds/ht_Prandtl

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%deltas(del)

#if NDIM == MDIM
   call Stencils_diffusion3d(solnData(HTN0_VAR, :, :, :), &
                             solnData(TEMP_VAR, :, :, :), &
                             del(DIR_X), del(DIR_Y), del(DIR_Z), &
                             diffusion_coeff, &
                             GRID_ILO, GRID_IHI, &
                             GRID_JLO, GRID_JHI, &
                             GRID_KLO, GRID_KHI)
#else
   call Stencils_diffusion2d(solnData(HTN0_VAR, :, :, :), &
                             solnData(TEMP_VAR, :, :, :), &
                             del(DIR_X), del(DIR_Y), &
                             diffusion_coeff, &
                             GRID_ILO, GRID_IHI, &
                             GRID_JLO, GRID_JHI)
#endif

   call tileDesc%releaseDataPtr(solnData, CENTER)

   call Timers_stop("HeatAD_diffusion")

end subroutine HeatAD_diffusion
