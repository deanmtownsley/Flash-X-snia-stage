!!***if* source/physics/sourceTerms/Outlet/OutletMain/phaseForcing/Outlet_setForcing
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

subroutine Outlet_setForcing(tileDesc, dt)

   use Outlet_data, ONLY: out_sink, out_flag, &
                          out_buffer, out_growthRate, &
                          out_QAuxLiq, out_QAuxGas, out_volAuxLiq, &
                          out_volAuxGas, out_QOutLiq, out_QOutGas, out_velRefScale, &
                          out_meshMe, out_xMin, out_xMax, out_yMin, out_yMax
#if NDIM == MDIM
   use Outlet_data, ONLY: out_zMin, out_zMax
#endif

   use Grid_interface, ONLY: Grid_getCellCoords
   use Grid_tile, ONLY: Grid_tile_t

   use out_interface, ONLY: out_lsDamping, out_velFrcPhased

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

   call Timers_start("Outlet_setForcing")

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
   call out_lsDamping(solnData(DFRC_VAR, :, :, :), &
                      solnData(DFUN_VAR, :, :, :), &
                      xCenter, yCenter, zCenter, boundBox, &
                      dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                      GRID_ILO, GRID_IHI, &
                      GRID_JLO, GRID_JHI, &
                      GRID_KLO, GRID_KHI, &
                      out_flag, out_sink, out_buffer, &
                      out_growthRate, &
                      out_xMin, out_xMax, out_yMin, out_yMax, 0., 0.)

   call out_velFrcPhased(facexData(VELC_FACE_VAR, :, :, :), &
                         facexData(VFRC_FACE_VAR, :, :, :), &
                         facexData(SIGM_FACE_VAR, :, :, :), &
                         solnData(DFUN_VAR, :, :, :), &
                         xCenter-del(IAXIS)/2, yCenter, zCenter, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         GRID_ILO, GRID_IHI+1, &
                         GRID_JLO, GRID_JHI, &
                         GRID_KLO, GRID_KHI, &
                         out_xMin, out_xMax, out_yMin, out_yMax, 0., 0., &
                         out_flag, out_buffer, out_growthRate, &
                         IAXIS, out_volAuxLiq, out_volAuxGas, out_QAuxLiq, out_QAuxGas, &
                         out_QOutLiq, out_QOutGas, out_velRefScale)

   call out_velFrcPhased(faceyData(VELC_FACE_VAR, :, :, :), &
                         faceyData(VFRC_FACE_VAR, :, :, :), &
                         faceyData(SIGM_FACE_VAR, :, :, :), &
                         solnData(DFUN_VAR, :, :, :), &
                         xCenter, yCenter-del(JAXIS)/2, zCenter, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         GRID_ILO, GRID_IHI, &
                         GRID_JLO, GRID_JHI+1, &
                         GRID_KLO, GRID_KHI, &
                         out_xMin, out_xMax, out_yMin, out_yMax, 0., 0., &
                         out_flag, out_buffer, out_growthRate, &
                         JAXIS, out_volAuxLiq, out_volAuxGas, out_QAuxLiq, out_QAuxGas, &
                         out_QOutLiq, out_QOutGas, out_velRefScale)

#else
   call tileDesc%getDataPtr(facezData, FACEZ)

   call out_lsDamping(solnData(DFRC_VAR, :, :, :), &
                      solnData(DFUN_VAR, :, :, :), &
                      xCenter, yCenter, zCenter, boundBox, &
                      dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                      GRID_ILO, GRID_IHI, &
                      GRID_JLO, GRID_JHI, &
                      GRID_KLO, GRID_KHI, &
                      out_flag, out_sink, out_buffer, &
                      out_growthRate, &
                      out_xMin, out_xMax, out_yMin, out_yMax, out_zMin, out_zMax)

   call out_velFrcPhased(facexData(VELC_FACE_VAR, :, :, :), &
                         facexData(VFRC_FACE_VAR, :, :, :), &
                         facexData(SIGM_FACE_VAR, :, :, :), &
                         solnData(DFUN_VAR, :, :, :), &
                         xCenter-del(IAXIS)/2, yCenter, zCenter, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         GRID_ILO, GRID_IHI+1, &
                         GRID_JLO, GRID_JHI, &
                         GRID_KLO, GRID_KHI, &
                         out_xMin, out_xMax, out_yMin, out_yMax, out_zMin, out_zMax, &
                         out_flag, out_buffer, out_growthRate, &
                         IAXIS, out_volAuxLiq, out_volAuxGas, out_QAuxLiq, out_QAuxGas, &
                         out_QOutLiq, out_QOutGas, out_velRefScale)

   call out_velFrcPhased(faceyData(VELC_FACE_VAR, :, :, :), &
                         faceyData(VFRC_FACE_VAR, :, :, :), &
                         faceyData(SIGM_FACE_VAR, :, :, :), &
                         solnData(DFUN_VAR, :, :, :), &
                         xCenter, yCenter-del(JAXIS)/2, zCenter, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         GRID_ILO, GRID_IHI, &
                         GRID_JLO, GRID_JHI+1, &
                         GRID_KLO, GRID_KHI, &
                         out_xMin, out_xMax, out_yMin, out_yMax, out_zMin, out_zMax, &
                         out_flag, out_buffer, out_growthRate, &
                         JAXIS, out_volAuxLiq, out_volAuxGas, out_QAuxLiq, out_QAuxGas, &
                         out_QOutLiq, out_QOutGas, out_velRefScale)

   call out_velFrcPhased(facezData(VELC_FACE_VAR, :, :, :), &
                         facezData(VFRC_FACE_VAR, :, :, :), &
                         facezData(SIGM_FACE_VAR, :, :, :), &
                         solnData(DFUN_VAR, :, :, :), &
                         xCenter, yCenter, zCenter-del(KAXIS)/2, &
                         dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                         GRID_ILO, GRID_IHI, &
                         GRID_JLO, GRID_JHI, &
                         GRID_KLO, GRID_KHI+1, &
                         out_xMin, out_xMax, out_yMin, out_yMax, out_zMin, out_zMax, &
                         out_flag, out_buffer, out_growthRate, &
                         KAXIS, out_volAuxLiq, out_volAuxGas, out_QAuxLiq, out_QAuxGas, &
                         out_QOutLiq, out_QOutGas, out_velRefScale)

   call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)

   call Timers_stop("Outlet_setForcing")

   return

end subroutine Outlet_setForcing
