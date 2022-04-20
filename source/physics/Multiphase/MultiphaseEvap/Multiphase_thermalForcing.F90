!!***if* source/physics/Multiphase/MultiphaseEvap/Multiphase_thermalForcing
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

subroutine Multiphase_thermalForcing(tileDesc)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t
   use mph_evapInterface, ONLY: mph_tempGfm2d, mph_tempGfm3d

!------------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:, :, :, :) :: solnData
   integer :: ierr, i, j, k
   real del(MDIM)
!------------------------------------------------------------------------------------------------
   nullify (solnData)

   call Timers_start("Multiphase_thermalForcing")

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%deltas(del)

#if NDIM < MDIM

   call mph_tempGfm2d(solnData(DFUN_VAR, :, :, :), &
                      solnData(NRMX_VAR, :, :, :), &
                      solnData(NRMY_VAR, :, :, :), &
                      (mph_invReynolds/mph_Prandtl)*solnData(mph_iAlphaCVar, :, :, :), &
                      solnData(mph_iTempVar, :, :, :), &
                      solnData(mph_iGfmVar, :, :, :), &
                      solnData(HFLQ_VAR, :, :, :), &
                      solnData(HFGS_VAR, :, :, :), &
                      mph_Tsat, &
                      del(DIR_X), del(DIR_Y), &
                      GRID_ILO_GC, GRID_IHI_GC, &
                      GRID_JLO_GC, GRID_JHI_GC)

#else

   call mph_tempGfm3d(solnData(DFUN_VAR, :, :, :), &
                      solnData(NRMX_VAR, :, :, :), &
                      solnData(NRMY_VAR, :, :, :), &
                      solnData(NRMZ_VAR, :, :, :), &
                      (mph_invReynolds/mph_Prandtl)*solnData(mph_iAlphaCVar, :, :, :), &
                      solnData(mph_iTempVar, :, :, :), &
                      solnData(mph_iGfmVar, :, :, :), &
                      solnData(HFLQ_VAR, :, :, :), &
                      solnData(HFGS_VAR, :, :, :), &
                      mph_Tsat, &
                      del(DIR_X), del(DIR_Y), del(DIR_Z), &
                      GRID_ILO_GC, GRID_IHI_GC, &
                      GRID_JLO_GC, GRID_JHI_GC, &
                      GRID_KLO_GC, GRID_KHI_GC)

#endif

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)

   call Timers_stop("Multiphase_thermalForcing")

   return

end subroutine Multiphase_thermalForcing
