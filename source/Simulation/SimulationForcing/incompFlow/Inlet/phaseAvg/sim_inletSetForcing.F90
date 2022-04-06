!!****f* source/Simulation/SimulationForcing/incompFlow/Inlet/phaseAvg/sim_inletSetForcing
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

subroutine sim_inletSetForcing(tileDesc, dt)

   use sim_outletData, ONLY: sim_outletFlag

   use sim_inletData, ONLY: sim_inletFlag, sim_inletGrowthRate, sim_inletBuffer, &
                            sim_inletSink

   use Simulation_data, ONLY: sim_meshMe, sim_xMin, sim_xMax, sim_yMin, sim_yMax
#if NDIM == MDIM
   use Simulation_data, ONLY: sim_zMin, sim_zMax
#endif

   use Grid_interface, ONLY: Grid_getCellCoords
   use Grid_tile, ONLY: Grid_tile_t

   use sim_inletInterface, ONLY: sim_inletLSDampingBlk2d, sim_inletLSDampingBlk3d, &
                                 sim_inletVelBlk2dPhased, sim_inletVelBlk3dPhased

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

   call Timers_start("sim_inletSetForcing")

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
   call sim_inletLSDampingBlk2d(solnData(DFRC_VAR, :, :, :), &
                                solnData(DFUN_VAR, :, :, :), &
                                xCenter, yCenter, boundBox, &
                                dt, del(IAXIS), del(JAXIS), &
                                GRID_ILO, GRID_IHI, &
                                GRID_JLO, GRID_JHI, &
                                sim_inletFlag, sim_inletSink, sim_inletBuffer, &
                                sim_inletGrowthRate, &
                                sim_xMin, sim_xMax, sim_yMin, sim_yMax)
#endif

   call sim_inletVelBlk2dPhased(facexData(VELC_FACE_VAR, :, :, :), &
                                faceyData(VELC_FACE_VAR, :, :, :), &
                                facexData(VFRC_FACE_VAR, :, :, :), &
                                faceyData(VFRC_FACE_VAR, :, :, :), &
                                solnData(DFUN_VAR, :, :, :), &
                                xCenter, yCenter, boundBox, &
                                dt, del(IAXIS), del(JAXIS), &
                                GRID_ILO, GRID_IHI, &
                                GRID_JLO, GRID_JHI, &
                                sim_inletFlag, sim_outletFlag, sim_inletBuffer, sim_inletGrowthRate, &
                                sim_xMin, sim_xMax, sim_yMin, sim_yMax, &
                                ins_gravX, ins_gravY)

#else
   call tileDesc%getDataPtr(facezData, FACEZ)

#ifdef MULTIPHASE_MAIN
   call sim_inletLSDampingBlk3d(solnData(DFRC_VAR, :, :, :), &
                                solnData(DFUN_VAR, :, :, :), &
                                xCenter, yCenter, zCenter, boundBox, &
                                dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                                GRID_ILO, GRID_IHI, &
                                GRID_JLO, GRID_JHI, &
                                GRID_KLO, GRID_KHI, &
                                sim_inletFlag, sim_inletSink, sim_inletBuffer, &
                                sim_inletGrowthRate, &
                                sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax)
#endif

   call sim_inletVelBlk3dPhased(facexData(VELC_FACE_VAR, :, :, :), &
                                faceyData(VELC_FACE_VAR, :, :, :), &
                                facezData(VELC_FACE_VAR, :, :, :), &
                                facexData(VFRC_FACE_VAR, :, :, :), &
                                faceyData(VFRC_FACE_VAR, :, :, :), &
                                facezData(VFRC_FACE_VAR, :, :, :), &
                                solnData(DFUN_VAR, :, :, :), &
                                xCenter, yCenter, zCenter, boundBox, &
                                dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                                GRID_ILO, GRID_IHI, &
                                GRID_JLO, GRID_JHI, &
                                GRID_KLO, GRID_KHI, &
                                sim_inletFlag, sim_outletFlag, sim_inletBuffer, sim_inletGrowthRate, &
                                sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax, &
                                ins_gravX, ins_gravY, ins_gravZ)

   call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)

   call Timers_stop("sim_inletSetForcing")

   return

end subroutine sim_inletSetForcing
