!!****if* source/Simulation/SimulationMain/incompFlow/CounterFlow/Simulation_initBlock
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
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID)
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified tile.
!!
!!  Reference:
!!
!!
!! ARGUMENTS
!!
!!  tile -          the tile to update
!!
!!
!!
!!
!!***
!!REORDER(4): solnData, face[xyz]Data

#include "constants.h"
#include "Simulation.h"

subroutine Simulation_initBlock(solnData, tileDesc)

   use Simulation_data
   use Grid_interface, ONLY: Grid_getCellCoords
   use Grid_tile, ONLY: Grid_tile_t

   implicit none

   !---Arguments ------------------------------------------------------------------------
   real, dimension(:, :, :, :), pointer :: solnData
   type(Grid_tile_t), intent(in)   :: tileDesc
   integer :: tileDescID

   !-------------------------------------------------------------------------------------
   integer, dimension(MDIM)       :: lo, hi
   real, allocatable, dimension(:) :: xGrid, yGrid, zGrid
   integer :: i, j, k
   real    :: xi, yi, zi
   real    :: del(MDIM)
   logical :: gcell = .true.
   real, pointer, dimension(:, :, :, :) :: facexData, faceyData, facezData
   real :: channelDepth
   real, parameter :: pi = acos(-1.0)

   !--------------------------------------------------------------------------------------
   nullify (facexData, faceyData, facezData)

   call tileDesc%deltas(del)
   lo = tileDesc%blkLimitsGC(LOW, :)
   hi = tileDesc%blkLimitsGC(HIGH, :)

   allocate (xGrid(lo(IAXIS):hi(IAXIS)))
   allocate (yGrid(lo(JAXIS):hi(JAXIS)))
   allocate (zGrid(lo(KAXIS):hi(KAXIS)))

   xGrid = 0.0
   yGrid = 0.0
   zGrid = 0.0

   call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xGrid)
   call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, lo, hi, yGrid)
#if NDIM == MDIM
   if (NDIM == 3) call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, lo, hi, zGrid)
#endif

   do k = lo(KAXIS), hi(KAXIS)
      do j = lo(JAXIS), hi(JAXIS)
         do i = lo(IAXIS), hi(IAXIS)
            xi = xGrid(i)
            yi = yGrid(j)
            zi = zGrid(k)

            !channelDepth = sim_channelDepth*(sim_yMax - yi)/sim_yMax
            !channelDepth = sim_channelDepth
            !
            !solnData(DFUN_VAR, i, j, k) = min(xi - (sim_xMin + channelDepth), (sim_xMax - channelDepth) - xi)

            solnData(DFUN_VAR, i, j, k) = min(xi-sim_xMin-sim_channelDepth-sim_nozzleAmp*cos(yi*pi/2), &
                                              sim_xMax-sim_channelDepth-xi+sim_nozzleAmp*cos(yi*pi/2+pi))

         end do
      end do
   end do

   deallocate (xGrid, yGrid, zGrid)

   allocate (xGrid(lo(IAXIS):hi(IAXIS)))
   allocate (yGrid(lo(JAXIS):hi(JAXIS)+1))
   allocate (zGrid(lo(KAXIS):hi(KAXIS)))

   xGrid = 0.0
   yGrid = 0.0
   zGrid = 0.0

   call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xGrid)
   call Grid_getCellCoords(JAXIS, FACES, tileDesc%level, lo, hi, yGrid)
#if NDIM == MDIM
   if (NDIM == 3) call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, lo, hi, zGrid)
#endif

   call tileDesc%getDataPtr(faceyData, FACEY)

   do k = lo(KAXIS), hi(KAXIS)
      do j = lo(JAXIS), hi(JAXIS)+1
         do i = lo(IAXIS), hi(IAXIS)

            yi = yGrid(j)

            if (0.5*(solnData(DFUN_VAR, i, j, k)+solnData(DFUN_VAR, i, j-1, k)) .gt. 0.0) then
               faceyData(VELC_FACE_VAR, i, j, k) = sim_gasFlowRate
            else
               faceyData(VELC_FACE_VAR, i, j, k) = sim_liqFlowRate
            end if

         end do
      end do
   end do

   call tileDesc%releaseDataPtr(faceyData, FACEY)

   deallocate (xGrid, yGrid, zGrid)

   return

end subroutine Simulation_initBlock
