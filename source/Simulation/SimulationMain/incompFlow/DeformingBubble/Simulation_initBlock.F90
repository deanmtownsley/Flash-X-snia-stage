!!****if* source/Simulation/SimulationMain/incompFlow/DeformingBubble/Simulation_initBlock
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

  !!$ Arguments -----------------------
   real, dimension(:, :, :, :), pointer :: solnData
   type(Grid_tile_t), intent(in) :: tileDesc
   integer :: tileDescID
  !!$ ---------------------------------

   integer :: i, j, k
   integer, dimension(MDIM) :: lo, hi
   real, allocatable, dimension(:) ::xCenter, yCenter, zCenter
   real :: xi, yi, zi
   logical :: gcell = .true.
   real, pointer, dimension(:, :, :, :) :: facexData, faceyData, facezData
   real, parameter :: pi = acos(-1.0)
   real :: del(MDIM)
   !----------------------------------------------------------------------
   nullify (facexData, faceyData, facezData)

   lo = tileDesc%blkLimitsGC(LOW, :)
   hi = tileDesc%blkLimitsGC(HIGH, :)

   allocate (xCenter(lo(IAXIS):hi(IAXIS)))
   allocate (yCenter(lo(JAXIS):hi(JAXIS)))
   allocate (zCenter(lo(KAXIS):hi(KAXIS)))

   xCenter = 0.0
   yCenter = 0.0
   zCenter = 0.0

   call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xCenter)
   if (NDIM >= 2) call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, lo, hi, yCenter)
   if (NDIM == 3) call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, lo, hi, zCenter)

   call tileDesc%deltas(del)
   call tileDesc%getDataPtr(facexData, FACEX)
   call tileDesc%getDataPtr(faceyData, FACEY)

   do k = lo(KAXIS), hi(KAXIS)
      do j = lo(JAXIS), hi(JAXIS)
         do i = lo(IAXIS), hi(IAXIS)
            xi = xCenter(i)
            yi = yCenter(j)

            solnData(DFUN_VAR, i, j, k) = 0.1 - sqrt((xi-0.75)**2 + (yi-0.75)**2)
         end do
      end do
   end do

   do k = lo(KAXIS), hi(KAXIS)
      do j = lo(JAXIS), hi(JAXIS)
         do i = lo(IAXIS), hi(IAXIS) + 1
            xi = xCenter(i) - del(IAXIS)/2
            yi = yCenter(j)

            facexData(VELC_FACE_VAR, i, j, k) = ((sin(pi*xi))**2)*sin(2*pi*yi)
         end do
      end do
   end do

   do k = lo(KAXIS), hi(KAXIS)
      do j = lo(JAXIS), hi(JAXIS) + 1
         do i = lo(IAXIS), hi(IAXIS)
            xi = xCenter(i)
            yi = yCenter(j) - del(JAXIS)/2

            faceyData(VELC_FACE_VAR, i, j, k) = -((sin(pi*yi))**2)*sin(2*pi*xi)
         end do
      end do
   end do


   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)
   deallocate (xCenter, yCenter, zCenter)

end subroutine Simulation_initBlock
