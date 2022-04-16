!!****f* source/Simulation/SimulationForcing/incompFlow/Outlet/sim_outletSetForcing
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

subroutine sim_outletSetForcing(tileDesc, dt)

   use sim_outletData, ONLY: sim_QOut, sim_outletSink, sim_outletFlag, &
                             sim_outletBuffer, sim_outletGrowthRate, sim_QAux, &
                             sim_volAux

   use Simulation_data, ONLY: sim_meshMe, sim_xMin, sim_xMax, sim_yMin, sim_yMax
#if NDIM == MDIM
   use Simulation_data, ONLY: sim_zMin, sim_zMax
#endif

   use Grid_interface, ONLY: Grid_getCellCoords
   use Grid_tile, ONLY: Grid_tile_t

   use sim_outletInterface, ONLY: sim_outletLSDamping, sim_outletVelFrc

   use IncompNS_data, ONLY: ins_gravX, ins_gravY, ins_gravZ
   use IncompNS_interface, ONLY: IncompNS_setVectorProp
   use Timers_interface, ONLY: Timers_start, Timers_stop

   implicit none
   include "Flashx_mpi.h"
   real, intent(in) :: dt
   type(Grid_tile_t), intent(in) :: tileDesc

!----------------------------------------------------------------------------------------
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   integer, dimension(2, MDIM)        :: blkLimits, blkLimitsGC
   integer, dimension(MDIM)          :: lo, hi
   real, dimension(GRID_IHI_GC)      :: xCenter
   real, dimension(GRID_JHI_GC)      :: yCenter
   real, dimension(GRID_KHI_GC)      :: zCenter
   real    :: del(MDIM)
   real    :: boundBox(LOW:HIGH, 1:MDIM)
   integer :: ierr

!----------------------------------------------------------------------------------------
   nullify (solnData, facexData, faceyData, facezData)

   call Timers_start("sim_outletSetForcing")

   blkLimits = tileDesc%limits
   blkLimitsGC = tileDesc%blkLimitsGC

   call tileDesc%deltas(del)
   call tileDesc%boundBox(boundBox)
   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%getDataPtr(facexData, FACEX)
   call tileDesc%getDataPtr(faceyData, FACEY)

   lo = blkLimitsGC(LOW, :)
   hi = blkLimitsGC(HIGH, :)

   xCenter = 0.0
   yCenter = 0.0
   zCenter = 0.0

   call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xCenter)
   call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, lo, hi, yCenter)

   if (NDIM == MDIM) call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, lo, hi, zCenter)

#if NDIM < MDIM

#ifdef MULTIPHASE_MAIN
   call sim_outletLSDamping(solnData(DFRC_VAR, :, :, :), &
                            solnData(DFUN_VAR, :, :, :), &
                            xCenter, yCenter, zCenter, boundBox, &
                            dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                            GRID_ILO, GRID_IHI, &
                            GRID_JLO, GRID_JHI, &
                            GRID_KLO, GRID_KHI, &
                            sim_outletFlag, sim_outletSink, sim_outletBuffer, &
                            sim_outletGrowthRate, &
                            sim_xMin, sim_xMax, sim_yMin, sim_yMax, 0., 0.)
#endif

   call sim_outletVelFrc(facexData(VELC_FACE_VAR, :, :, :), &
                         facexData(VFRC_FACE_VAR, :, :, :), &
                         xCenter - del(IAXIS)/2, yCenter, zCenter, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         GRID_ILO, GRID_IHI + 1, &
                         GRID_JLO, GRID_JHI, &
                         GRID_KLO, GRID_KHI, &
                         sim_xMin, sim_xMax, sim_yMin, sim_yMax, 0., 0., &
                         sim_outletFlag, sim_outletBuffer, sim_outletGrowthRate, &
                         IAXIS, sim_volAux, sim_QAux, sim_QOut)

   call sim_outletVelFrc(faceyData(VELC_FACE_VAR, :, :, :), &
                         faceyData(VFRC_FACE_VAR, :, :, :), &
                         xCenter, yCenter - del(JAXIS)/2, zCenter, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         GRID_ILO, GRID_IHI, &
                         GRID_JLO, GRID_JHI + 1, &
                         GRID_KLO, GRID_KHI, &
                         sim_xMin, sim_xMax, sim_yMin, sim_yMax, 0., 0., &
                         sim_outletFlag, sim_outletBuffer, sim_outletGrowthRate, &
                         JAXIS, sim_volAux, sim_QAux, sim_QOut)

#else
   call tileDesc%getDataPtr(facezData, FACEZ)

#ifdef MULTIPHASE_MAIN
   call sim_outletLSDamping(solnData(DFRC_VAR, :, :, :), &
                            solnData(DFUN_VAR, :, :, :), &
                            xCenter, yCenter, zCenter, boundBox, &
                            dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                            GRID_ILO, GRID_IHI, &
                            GRID_JLO, GRID_JHI, &
                            GRID_KLO, GRID_KHI, &
                            sim_outletFlag, sim_outletSink, sim_outletBuffer, &
                            sim_outletGrowthRate, &
                            sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax)
#endif

   call sim_outletVelFrc(facexData(VELC_FACE_VAR, :, :, :), &
                         facexData(VFRC_FACE_VAR, :, :, :), &
                         xCenter - del(IAXIS)/2, yCenter, zCenter, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         GRID_ILO, GRID_IHI + 1, &
                         GRID_JLO, GRID_JHI, &
                         GRID_KLO, GRID_KHI, &
                         sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax, &
                         sim_outletFlag, sim_outletBuffer, sim_outletGrowthRate, &
                         IAXIS, sim_volAux, sim_QAux, sim_QOut)

   call sim_outletVelFrc(faceyData(VELC_FACE_VAR, :, :, :), &
                         faceyData(VFRC_FACE_VAR, :, :, :), &
                         xCenter, yCenter - del(JAXIS)/2, zCenter, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         GRID_ILO, GRID_IHI, &
                         GRID_JLO, GRID_JHI + 1, &
                         GRID_KLO, GRID_KHI, &
                         sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax, &
                         sim_outletFlag, sim_outletBuffer, sim_outletGrowthRate, &
                         JAXIS, sim_volAux, sim_QAux, sim_QOut)

   call sim_outletVelFrc(facezData(VELC_FACE_VAR, :, :, :), &
                         facezData(VFRC_FACE_VAR, :, :, :), &
                         xCenter, yCenter, zCenter - del(KAXIS)/2, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         GRID_ILO, GRID_IHI, &
                         GRID_JLO, GRID_JHI, &
                         GRID_KLO, GRID_KHI + 1, &
                         sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax, &
                         sim_outletFlag, sim_outletBuffer, sim_outletGrowthRate, &
                         KAXIS, sim_volAux, sim_QAux, sim_QOut)

   call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)

   call Timers_stop("sim_outletSetForcing")

   return

end subroutine sim_outletSetForcing
