!!***if* source/physics/ImBound/ImBoundMain/ImBound_mapToGrid
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
#include "Simulation.h"

subroutine ImBound_mapToGrid(tileDesc, bodyInfo)

   use ImBound_data
   use ImBound_type, ONLY: ImBound_type_t
   use ib_interface, ONLY: ib_bruteForceMap, ib_annMap
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t
   use Grid_interface, ONLY: Grid_getCellCoords
   use Stencils_interface, ONLY: Stencils_lsNormals2d, Stencils_lsNormals3d, &
                                 Stencils_lsCurvature2d, Stencils_lsCurvature3d

!-----------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc
   type(ImBound_type_t), intent(in) :: bodyInfo

   real, pointer, dimension(:, :, :, :) :: solnData
   real :: del(MDIM)
   real, dimension(GRID_IHI_GC) :: xCenter
   real, dimension(GRID_JHI_GC) :: yCenter
   real, dimension(GRID_KHI_GC) :: zCenter
   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   integer, dimension(MDIM) :: lo, hi
   real, dimension(LOW:HIGH, 1:MDIM) :: boundBox

!-----------------------------------------------------------------------------------------
   nullify (solnData)

   call Timers_start("ImBound_mapToGrid")

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%deltas(del)
   call tileDesc%boundBox(boundBox)

   blkLimits = tileDesc%limits
   blkLimitsGC = tileDesc%blkLimitsGC

   lo = blkLimitsGC(LOW, :)
   hi = blkLimitsGC(HIGH, :)

   xCenter = 0.0
   yCenter = 0.0
   zCenter = 0.0

   call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xCenter)
   call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, lo, hi, yCenter)

   if (NDIM == MDIM) call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, lo, hi, zCenter)

#if NDIM == MDIM
   if (ib_bruteForceMapping) then
      call ib_bruteForceMap(solnData(LMDA_VAR, :, :, :), &
                            xCenter, yCenter, zCenter, &
                            del(IAXIS), del(JAXIS), del(KAXIS), &
                            GRID_ILO_GC, GRID_IHI_GC, &
                            GRID_JLO_GC, GRID_JHI_GC, &
                            GRID_KLO_GC, GRID_KHI_GC, &
                            bodyInfo)

   else
      call ib_annMap(solnData(LMDA_VAR, :, :, :), &
                     xCenter, yCenter, zCenter, &
                     del(IAXIS), del(JAXIS), del(KAXIS), &
                     GRID_ILO_GC, GRID_IHI_GC, &
                     GRID_JLO_GC, GRID_JHI_GC, &
                     GRID_KLO_GC, GRID_KHI_GC, &
                     bodyInfo)

   end if

   call Stencils_lsNormals3d(solnData(LMDA_VAR, :, :, :), &
                             solnData(NMLX_VAR, :, :, :), &
                             solnData(NMLY_VAR, :, :, :), &
                             solnData(NMLZ_VAR, :, :, :), &
                             del(IAXIS), del(JAXIS), del(KAXIS), &
                             GRID_ILO_GC, GRID_IHI_GC, &
                             GRID_JLO_GC, GRID_JHI_GC, &
                             GRID_KLO_GC, GRID_KHI_GC)

   call Stencils_lsCurvature3d(solnData(LCRV_VAR, :, :, :), &
                               solnData(LMDA_VAR, :, :, :), &
                               solnData(NMLX_VAR, :, :, :), &
                               solnData(NMLY_VAR, :, :, :), &
                               solnData(NMLZ_VAR, :, :, :), &
                               del(IAXIS), del(JAXIS), del(KAXIS), &
                               GRID_ILO_GC, GRID_IHI_GC, &
                               GRID_JLO_GC, GRID_JHI_GC, &
                               GRID_KLO_GC, GRID_KHI_GC)

#else
   if (ib_bruteForceMapping) then
      call ib_bruteForceMap(solnData(LMDA_VAR, :, :, :), &
                            xCenter, yCenter, &
                            del(IAXIS), del(JAXIS), &
                            GRID_ILO_GC, GRID_IHI_GC, &
                            GRID_JLO_GC, GRID_JHI_GC, &
                            bodyInfo)

   else
      call ib_annMap(solnData(LMDA_VAR, :, :, :), &
                     xCenter, yCenter, &
                     del(IAXIS), del(JAXIS), &
                     GRID_ILO_GC, GRID_IHI_GC, &
                     GRID_JLO_GC, GRID_JHI_GC, &
                     bodyInfo)

   end if

   call Stencils_lsNormals2d(solnData(LMDA_VAR, :, :, :), &
                             solnData(NMLX_VAR, :, :, :), &
                             solnData(NMLY_VAR, :, :, :), &
                             del(IAXIS), del(JAXIS), &
                             GRID_ILO_GC, GRID_IHI_GC, &
                             GRID_JLO_GC, GRID_JHI_GC)

   call Stencils_lsCurvature2d(solnData(LCRV_VAR, :, :, :), &
                               solnData(LMDA_VAR, :, :, :), &
                               solnData(NMLX_VAR, :, :, :), &
                               solnData(NMLY_VAR, :, :, :), &
                               del(IAXIS), del(JAXIS), &
                               GRID_ILO_GC, GRID_IHI_GC, &
                               GRID_JLO_GC, GRID_JHI_GC)
#endif

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call Timers_stop("ImBound_mapToGrid")

end subroutine ImBound_mapToGrid
