!!***if* source/physics/Multiphase/MultiphaseMain/Multiphase_setFluidProps
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
!!REORDER(4): solnData, face[xyz]Data

#include "Simulation.h"
#include "constants.h"
#include "Multiphase.h"

subroutine Multiphase_setFluidProps(tileDesc)

   use Multiphase_data
   use Stencils_interface, ONLY: Stencils_lsCenterProps, &
                                 Stencils_lsFaceProps2d, &
                                 Stencils_lsFaceProps3d
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t

!---------------------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   integer, dimension(2, MDIM) :: stnLimits = 1, stnLimitsGC = 1
   logical :: gcMask(NUNK_VARS + NDIM*NFACE_VARS)
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   integer :: ierr, i, j, k
   real del(MDIM)
   real minCellDiag

!----------------------------------------------------------------------------------------------------------
   nullify (solnData, facexData, faceyData, facezData)

   call Timers_start("Multiphase_setFluidProps")

   stnLimitsGC(LOW, 1:NDIM) = tileDesc%limits(LOW, 1:NDIM) - tileDesc%blkLimitsGC(LOW, 1:NDIM) + 1 - NGUARD
   stnLimitsGC(HIGH, 1:NDIM) = tileDesc%limits(HIGH, 1:NDIM) - tileDesc%blkLimitsGC(LOW, 1:NDIM) + 1 + NGUARD

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%getDataPtr(facexData, FACEX)
   call tileDesc%getDataPtr(faceyData, FACEY)
#if NDIM == MDIM
   call tileDesc%getDataPtr(facezData, FACEZ)
#endif
   call tileDesc%deltas(del)

   minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.+del(DIR_Z)**2.)

#if NDIM < MDIM
   call Stencils_lsFaceProps2d(solnData(DFUN_VAR, :, :, :), &
                               facexData(mph_iRhoFVar, :, :, :), &
                               faceyData(mph_iRhoFVar, :, :, :), &
                               1./mph_rhoGas, &
                               stnLimitsGC(LOW, IAXIS), stnLimitsGC(HIGH, IAXIS), &
                               stnLimitsGC(LOW, JAXIS), stnLimitsGC(HIGH, JAXIS), &
                               iSmear=mph_iPropSmear*minCellDiag)

#else
   call Stencils_lsFaceProps3d(solnData(DFUN_VAR, :, :, :), &
                               facexData(mph_iRhoFVar, :, :, :), &
                               faceyData(mph_iRhoFVar, :, :, :), &
                               facezData(mph_iRhoFVar, :, :, :), &
                               1./mph_rhoGas, &
                               stnLimitsGC(LOW, IAXIS), stnLimitsGC(HIGH, IAXIS), &
                               stnLimitsGC(LOW, JAXIS), stnLimitsGC(HIGH, JAXIS), &
                               stnLimitsGC(LOW, KAXIS), stnLimitsGC(HIGH, KAXIS), &
                               iSmear=mph_iPropSmear*minCellDiag)
#endif

   call Stencils_lsCenterProps(solnData(DFUN_VAR, :, :, :), &
                               solnData(mph_iMuCVar, :, :, :), &
                               mph_muGas, &
                               stnLimitsGC(LOW, IAXIS), stnLimitsGC(HIGH, IAXIS), &
                               stnLimitsGC(LOW, JAXIS), stnLimitsGC(HIGH, JAXIS), &
                               stnLimitsGC(LOW, KAXIS), stnLimitsGC(HIGH, KAXIS), &
                               iSmear=mph_iPropSmear*minCellDiag)

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)
#if NDIM == MDIM
   call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif
   call Timers_stop("Multiphase_setFluidProps")

   return
end subroutine Multiphase_setFluidProps
