!!***if* source/physics/Multiphase/MultiphaseMain/Multiphase_indicators
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

#include "Simulation.h"
#include "constants.h"
#include "Multiphase.h"

subroutine Multiphase_indicators()

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_interface, ONLY: Grid_getTileIterator, Grid_releaseTileIterator
   use Grid_tile, ONLY: Grid_tile_t
   use Grid_iterator, ONLY: Grid_iterator_t
   use Stencils_interface, ONLY: Stencils_lsRedistance2d, Stencils_lsRedistance3d

!-----------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:, :, :, :) :: solnData
   integer :: ierr, i, j, k
   real del(MDIM)
   real :: volSum, volSumAll
   type(Grid_tile_t) :: tileDesc
   type(Grid_iterator_t) :: itor

!-----------------------------------------------------------------------------------------
   nullify (solnData)

   volSum = 0.0
   volSumAll = 0.0

   call Grid_getTileIterator(itor, nodetype=LEAF)
   do while (itor%isValid())
      call itor%currentTile(tileDesc)
      call tileDesc%getDataPtr(solnData, CENTER)
      blkLimits = tileDesc%limits
      call tileDesc%deltas(del)
      do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
         do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
            do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
               if (solnData(DFUN_VAR, i, j, k) .gt. 0) then
#if NDIM < MDIM
                  volSum = volSum + del(DIR_X)*del(DIR_Y)
#else
                  volSum = volSum + del(DIR_X)*del(DIR_Y)*del(DIR_Z)
#endif
               end if
            end do
         end do
      end do
      call tileDesc%releaseDataPtr(solnData, CENTER)
      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)

   call MPI_Allreduce(volSum, volSumAll, 1, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)
   if (mph_meshMe .eq. 0) print *, "----------------------------------------"
   if (mph_meshMe .eq. 0) print *, "Total Gas Volume: ", volSumAll
   if (mph_meshMe .eq. 0) print *, "----------------------------------------"

   return

end subroutine Multiphase_indicators
