!!****f* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterCheckSites
!!
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
!!***
!!REORDER(4): solnData, face[xyz]Data

#include "Simulation.h"
#include "constants.h"

subroutine sim_heaterCheckSites(tileDesc)

   use Grid_interface, ONLY: Grid_getCellCoords
   use Grid_tile, ONLY: Grid_tile_t
   use sim_heaterInterface, ONLY: sim_heaterCheckSitesBlk2d, sim_heaterCheckSitesBlk3d
   use Timers_interface, ONLY: Timers_start, Timers_stop

!----------------------------------------------------------------------------------------
   implicit none
   type(Grid_tile_t), intent(in) :: tileDesc

   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   integer, dimension(2, MDIM)        :: blkLimits, blkLimitsGC
   integer, dimension(MDIM)          :: lo, hi
   real, dimension(GRID_IHI_GC)      :: xCenter
   real, dimension(GRID_JHI_GC)      :: yCenter
   real, dimension(GRID_KHI_GC)      :: zCenter
   real    :: del(MDIM)
   real    :: boundBox(LOW:HIGH, 1:MDIM)

!----------------------------------------------------------------------------------------
   nullify (solnData, facexData, faceyData, facezData)

   call Timers_start("sim_heaterCheckSites")

#ifdef MULTIPHASE_EVAPORATION
   blkLimits = tileDesc%limits
   blkLimitsGC = tileDesc%blkLimitsGC
   call tileDesc%deltas(del)
   call tileDesc%boundBox(boundBox)
   call tileDesc%getDataPtr(solnData, CENTER)

   lo = blkLimitsGC(LOW, :)
   hi = blkLimitsGC(HIGH, :)

   xCenter = 0.0
   yCenter = 0.0
   zCenter = 0.0
   call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xCenter)
   call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, lo, hi, yCenter)
   if (NDIM == MDIM) call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, lo, hi, zCenter)

#if NDIM < MDIM
   call sim_heaterCheckSitesBlk2d(solnData(DFUN_VAR, :, :, :), &
                                  xCenter, yCenter, boundBox, &
                                  GRID_ILO_GC, GRID_IHI_GC, &
                                  GRID_JLO_GC, GRID_JHI_GC)
#else
   call sim_heaterCheckSitesBlk3d(solnData(DFUN_VAR, :, :, :), &
                                  xCenter, yCenter, zCenter, boundBox, &
                                  GRID_ILO_GC, GRID_IHI_GC, &
                                  GRID_JLO_GC, GRID_JHI_GC, &
                                  GRID_KLO_GC, GRID_KHI_GC)
#endif

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
#endif

   call Timers_stop("sim_heaterCheckSites")

   return
end subroutine sim_heaterCheckSites
