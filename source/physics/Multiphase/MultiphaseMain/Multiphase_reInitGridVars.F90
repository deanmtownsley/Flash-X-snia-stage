!!****if* source/physics/Multiphase/MultiphaseMain/Multiphase_reInitGridVars
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
!!
!! Multiphase_reInitGridVars
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!***
!!REORDER(4): solnData

#include "constants.h"
#include "Multiphase.h"
#include "Simulation.h"

subroutine Multiphase_reInitGridVars(solnData, logc, higc)

   use Timers_interface, ONLY: Timers_start, Timers_stop

   !------------------------------------------------------------------------------------------
   implicit none
   integer, dimension(MDIM),intent(IN) :: logc, higc
   real, pointer, dimension(:, :, :, :) :: solnData
   integer :: i, j, k
   !------------------------------------------------------------------------------------------

   call Timers_start("Multiphase_reInitGridVars")


!!$   do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
!!$      do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
!!$         do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)

            ! DEVNOTE (10/24/2023):
            ! Commenting intialization see
            ! explanation below
            !solnData(CURV_VAR, i, j, k) = 0

   solnData(DFRC_VAR, :,:,:) = 0.
   solnData(HDN0_VAR, :,:,:) = 0.
   
            ! DEVNOTE (10/24/2023):
            ! Commenting initializations that are
            ! unnecessary to improve performance
            ! during high blocks/process loading
#ifdef MULTIPHASE_EVAPORATION
            !solnData(HFLQ_VAR, :,:,:) = 0.
            !solnData(HFGS_VAR, :,:,:) = 0.
#endif

            !solnData(NRMX_VAR, :,:,:) = 0.
            !solnData(NRMY_VAR, :,:,:) = 0.

#if NDIM == MDIM
            !solnData(NRMZ_VAR, :,:,:) = 0.
#endif

   solnData(PFUN_VAR, :,:,:) = 0.
   solnData(SMHV_VAR, :,:,:) = 0.
            
!!$         end do
!!$      end do
!!$   end do

   ! Release pointers:

   call Timers_stop("Multiphase_reInitGridVars")

   return
end subroutine Multiphase_reInitGridVars
