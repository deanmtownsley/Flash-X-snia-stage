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

subroutine Multiphase_redistance(solnData,del,lo,hi, iteration)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Stencils_interface, ONLY: Stencils_lsRedistance2d, Stencils_lsRedistance3d

!-----------------------------------------------------------------------------------------
   implicit none
   real, pointer, dimension(:, :, :, :) :: solnData
   integer, intent(in) :: iteration
   integer,dimension(MDIM), intent(IN) :: lo, hi
   real,dimension(MDIM),intent(IN) :: del

   real :: lsDT, minCellDiag
!-----------------------------------------------------------------------------------------

   call Timers_start("Multiphase_redistance")

   if (iteration .eq. 1) then
      solnData(HDN0_VAR, :, :, :) = solnData(DFUN_VAR, :, :, :)
   end if

   minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.+del(DIR_Z)**2)
   lsDT = minCellDiag/5.0d0
#if NDIM < MDIM
   !--------------------------------------------
   ! Call DFUN re-initialization routine for 2D:
   !--------------------------------------------
   call Stencils_lsRedistance2d(solnData(DFUN_VAR, :, :, :), &
                                solnData(HDN0_VAR, :, :, :), &
                                lsDT, del(DIR_X), del(DIR_Y), &
                                lo( IAXIS), hi( IAXIS), &
                                lo( JAXIS), hi( JAXIS))
#else
   call Stencils_lsRedistance3d(solnData(DFUN_VAR, :, :, :), &
                                solnData(HDN0_VAR, :, :, :), &
                                lsDT, del(DIR_X), del(DIR_Y), del(DIR_Z), &
                                lo( IAXIS), hi( IAXIS), &
                                lo( JAXIS), hi( JAXIS), &
                                lo( KAXIS), hi( KAXIS))
#endif
   call Timers_stop("Multiphase_redistance")

   return

end subroutine Multiphase_redistance
