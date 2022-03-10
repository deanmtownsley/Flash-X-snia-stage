!!***if* source/physics/Multiphase/MultiphaseEvap/Multiphase_divergence
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
!!REORDER(4): solnData,face[xyz]Data

#include "constants.h"
#include "Multiphase.h"
#include "Simulation.h"

subroutine Multiphase_divergence(tileDesc)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t
   use Stencils_interface, ONLY: Stencils_cnt_advectUpwind2d, Stencils_cnt_advectUpwind3d
   use mph_evapInterface, ONLY: mph_evapDivergence2d, mph_evapDivergence3d

!------------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   integer :: ierr, i, j, k
   real del(MDIM)
!------------------------------------------------------------------------------------------------
   nullify (solnData, facexData, faceyData, facezData)

   call Timers_start("Multiphase_divergence")

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%getDataPtr(facexData, FACEX)
   call tileDesc%getDataPtr(faceyData, FACEY)
   call tileDesc%deltas(del)

#if NDIM < MDIM
   call mph_evapDivergence2d(solnData(mph_iDivCvar, :, :, :), &
                             facexData(mph_iRhoFVar, :, :, :), &
                             faceyData(mph_iRhoFVar, :, :, :), &
                             solnData(NRMX_VAR, :, :, :), &
                             solnData(NRMY_VAR, :, :, :), &
                             solnData(MFLX_VAR, :, :, :), &
                             del(DIR_X), del(DIR_Y), &
                             GRID_ILO, GRID_IHI, &
                             GRID_JLO, GRID_JHI)

#else

   call tileDesc%getDataPtr(facezData, FACEZ)

   call mph_evapDivergence3d(solnData(mph_iDivCvar, :, :, :), &
                             facexData(mph_iRhoFVar, :, :, :), &
                             faceyData(mph_iRhoFVar, :, :, :), &
                             facezData(mph_iRhoFVar, :, :, :), &
                             solnData(NRMX_VAR, :, :, :), &
                             solnData(NRMY_VAR, :, :, :), &
                             solnData(NRMZ_VAR, :, :, :), &
                             solnData(MFLX_VAR, :, :, :), &
                             del(DIR_X), del(DIR_Y), del(DIR_Z), &
                             GRID_ILO, GRID_IHI, &
                             GRID_JLO, GRID_JHI, &
                             GRID_KLO, GRID_KHI)

   call tileDesc%releaseDataPtr(facezData, FACEZ)

#endif

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)

   call Timers_stop("Multiphase_divergence")

   return

end subroutine Multiphase_divergence
