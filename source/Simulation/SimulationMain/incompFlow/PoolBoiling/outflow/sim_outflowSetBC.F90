!!****f* source/Simulation/SimulationMain/incompFlow/PoolBoiling/localAPI/outflow/sim_outflowSetBC
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

subroutine sim_outflowSetBC(dt)

   use sim_outflowData, ONLY: sim_outflowVel
   use Simulation_data, ONLY: sim_meshMe
   use Grid_interface, ONLY: Grid_getTileIterator, Grid_releaseTileIterator, Grid_getCellCoords
   use Grid_tile, ONLY: Grid_tile_t
   use Grid_iterator, ONLY: Grid_iterator_t
   use sim_outflowInterface, ONLY: sim_outflowLSDampingBlk2d, sim_outflowLSDampingBlk3d, &
                                   sim_outflowVelBlk2d, sim_outflowVelBlk3d
   use IncompNS_interface, ONLY: IncompNS_setVectorProp
   use Timers_interface, ONLY: Timers_start, Timers_stop

   implicit none
   include "Flashx_mpi.h"
   real, intent(in) :: dt

!----------------------------------------------------------------------------------------
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   integer, dimension(2, MDIM)        :: blkLimits, blkLimitsGC
   integer, dimension(MDIM)          :: lo, hi
   real, dimension(GRID_IHI_GC)      :: xCenter
   real, dimension(GRID_JHI_GC)      :: yCenter
   real, dimension(GRID_KHI_GC)      :: zCenter
   real    :: del(MDIM)
   type(Grid_tile_t) :: tileDesc
   type(Grid_iterator_t) :: itor
   real    :: boundBox(LOW:HIGH, 1:MDIM)
   real    :: velOut
   integer :: ierr

!----------------------------------------------------------------------------------------
   nullify (solnData, facexData, faceyData, facezData)

   call Timers_start("sim_outflowSetBC")

   velOut = 1.

   call Grid_getTileIterator(itor, nodetype=LEAF)
   !
   do while (itor%isValid())
      call itor%currentTile(tileDesc)
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

      call sim_outflowLSDampingBlk2d(solnData(DFRC_VAR, :, :, :), &
                                     solnData(DFUN_VAR, :, :, :), &
                                     xCenter, yCenter, boundBox, &
                                     dt, del(IAXIS), del(JAXIS), &
                                     GRID_ILO_GC, GRID_IHI_GC, &
                                     GRID_JLO_GC, GRID_JHI_GC)

      call sim_outflowVelBlk2d(velOut, &
                               facexData(VELC_FACE_VAR, :, :, :), &
                               faceyData(VELC_FACE_VAR, :, :, :), &
                               facexData(VFRC_FACE_VAR, :, :, :), &
                               faceyData(VFRC_FACE_VAR, :, :, :), &
                               xCenter, yCenter, boundBox, &
                               dt, del(IAXIS), del(JAXIS), &
                               GRID_ILO_GC, GRID_IHI_GC, &
                               GRID_JLO_GC, GRID_JHI_GC)

#if NDIM == MDIM
      call tileDesc%getDataPtr(facezData, FACEZ)

      call sim_outflowLSDampingBlk3d(solnData(DFRC_VAR, :, :, :), &
                                     solnData(DFUN_VAR, :, :, :), &
                                     xCenter, yCenter, zCenter, boundBox, &
                                     dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                                     GRID_ILO_GC, GRID_IHI_GC, &
                                     GRID_JLO_GC, GRID_JHI_GC, &
                                     GRID_KLO_GC, GRID_KHI_GC)

      call sim_outflowVelBlk3d(velOut, &
                               facexData(VELC_FACE_VAR, :, :, :), &
                               faceyData(VELC_FACE_VAR, :, :, :), &
                               facezData(VELC_FACE_VAR, :, :, :), &
                               facexData(VFRC_FACE_VAR, :, :, :), &
                               faceyData(VFRC_FACE_VAR, :, :, :), &
                               facezData(VFRC_FACE_VAR, :, :, :), &
                               xCenter, yCenter, zCenter, boundBox, &
                               dt, del(IAXIS), del(JAXIS), del(KAXIS), &
                               GRID_ILO_GC, GRID_IHI_GC, &
                               GRID_JLO_GC, GRID_JHI_GC, &
                               GRID_KLO_GC, GRID_KHI_GC)

      call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

      ! Release pointers:
      call tileDesc%releaseDataPtr(solnData, CENTER)
      call tileDesc%releaseDataPtr(facexData, FACEX)
      call tileDesc%releaseDataPtr(faceyData, FACEY)
      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)

   sim_outflowVel = 0.

   call MPI_Allreduce(velOut, sim_outflowVel(HIGH, JAXIS), 1, FLASH_REAL, &
                      MPI_MAX, MPI_COMM_WORLD, ierr)

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'Outflow Velocity Low  =', sim_outflowVel(LOW, :)
      write (*, *) 'Outflow Velocity High =', sim_outflowVel(HIGH, :)
   end if

   call IncompNS_setVectorProp("Outflow_Vel_Low", sim_outflowVel(LOW, :))
   call IncompNS_setVectorProp("Outflow_Vel_High", sim_outflowVel(HIGH, :))

   call Timers_stop("sim_outflowSetBC")

   return

end subroutine sim_outflowSetBC
