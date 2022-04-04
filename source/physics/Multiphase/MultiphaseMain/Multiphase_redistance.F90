!!***if* source/physics/Multiphase/MultiphaseMain/Multiphase_redistance
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

subroutine Multiphase_redistance(tileDesc, iteration)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t
   use Stencils_interface, ONLY: Stencils_lsRedistance2d, Stencils_lsRedistance3d

!-----------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   integer, intent(in) :: iteration
   type(Grid_tile_t), intent(in) :: tileDesc

   integer, dimension(2, MDIM) :: stnLimits = 1
   real, pointer, dimension(:, :, :, :) :: solnData
   integer :: ierr
   real del(MDIM)
   real :: lsDT, minCellDiag
!-----------------------------------------------------------------------------------------
   nullify (solnData)

   call Timers_start("Multiphase_redistance")

   call tileDesc%getDataPtr(solnData, CENTER)

   stnLimits(LOW, 1:NDIM) = tileDesc%limits(LOW, 1:NDIM) - tileDesc%blkLimitsGC(LOW, 1:NDIM) + 1
   stnLimits(HIGH, 1:NDIM) = tileDesc%limits(HIGH, 1:NDIM) - tileDesc%blkLimitsGC(LOW, 1:NDIM) + 1

   if (iteration .eq. 1) then
      solnData(HDN0_VAR, :, :, :) = solnData(DFUN_VAR, :, :, :)
   end if

   call tileDesc%deltas(del)
   minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.+del(DIR_Z)**2)
   lsDT = minCellDiag/5.0d0
#if NDIM < MDIM
   !--------------------------------------------
   ! Call DFUN re-initialization routine for 2D:
   !--------------------------------------------
   call Stencils_lsRedistance2d(solnData(DFUN_VAR, :, :, :), &
                                solnData(HDN0_VAR, :, :, :), &
                                lsDT, del(DIR_X), del(DIR_Y), &
                                stnLimits(LOW, IAXIS), stnLimits(HIGH, IAXIS), &
                                stnLimits(LOW, JAXIS), stnLimits(HIGH, JAXIS))
#else
   call Stencils_lsRedistance3d(solnData(DFUN_VAR, :, :, :), &
                                solnData(HDN0_VAR, :, :, :), &
                                lsDT, del(DIR_X), del(DIR_Y), del(DIR_Z), &
                                stnLimits(LOW, IAXIS), stnLimits(HIGH, IAXIS), &
                                stnLimits(LOW, JAXIS), stnLimits(HIGH, JAXIS), &
                                stnLimits(LOW, KAXIS), stnLimits(HIGH, KAXIS))
#endif
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call Timers_stop("Multiphase_redistance")

   return

end subroutine Multiphase_redistance
