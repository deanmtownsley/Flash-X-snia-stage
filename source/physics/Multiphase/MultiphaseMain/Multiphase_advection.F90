!!***if* source/physics/Multiphase/MultiphaseMain/Multiphase_advection
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

#include "constants.h"
#include "Multiphase.h"
#include "Simulation.h"

subroutine Multiphase_advection(tileDesc)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t
   use Stencils_interface, ONLY: Stencils_advectWeno2d, Stencils_advectWeno3d
   use mph_evapInterface, ONLY: mph_evapVelocity2d, mph_evapVelocity3d

!-----------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   integer, dimension(2, MDIM) :: stnLimits = 1
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   real del(MDIM)
!-----------------------------------------------------------------------------------------
   nullify (solnData, facexData, faceyData, facezData)

   call Timers_start("Multiphase_advection")

   stnLimits(LOW, 1:NDIM) = tileDesc%limits(LOW, 1:NDIM) - tileDesc%blkLimitsGC(LOW, 1:NDIM) + 1
   stnLimits(HIGH, 1:NDIM) = tileDesc%limits(HIGH, 1:NDIM) - tileDesc%blkLimitsGC(LOW, 1:NDIM) + 1

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%getDataPtr(facexData, FACEX)
   call tileDesc%getDataPtr(faceyData, FACEY)
#if NDIM == MDIM
   call tileDesc%getDataPtr(facezData, FACEZ)
#endif
   call tileDesc%deltas(del)

#if NDIM < MDIM

#ifdef MULTIPHASE_EVAPORATION
   call mph_evapVelocity2d(facexData(mph_iVelFVar, :, :, :), &
                           faceyData(mph_iVelFVar, :, :, :), &
                           facexData(mph_iRhoFVar, :, :, :), &
                           faceyData(mph_iRhoFVar, :, :, :), &
                           solnData(NRMX_VAR, :, :, :), &
                           solnData(NRMY_VAR, :, :, :), &
                           solnData(MFLX_VAR, :, :, :), &
                           stnLimits(LOW, IAXIS), stnLimits(HIGH, IAXIS), &
                           stnLimits(LOW, JAXIS), stnLimits(HIGH, JAXIS))
#endif

   call Stencils_advectWeno2d(solnData(HDN0_VAR, :, :, :), &
                              solnData(DFUN_VAR, :, :, :), &
                              facexData(mph_iVelFVar, :, :, :), &
                              faceyData(mph_iVelFVar, :, :, :), &
                              del(DIR_X), del(DIR_Y), &
                              stnLimits(LOW, IAXIS), stnLimits(HIGH, IAXIS), &
                              stnLimits(LOW, JAXIS), stnLimits(HIGH, JAXIS), &
                              center=1, facex=0, facey=0)

#ifdef MULTIPHASE_EVAPORATION
   call mph_evapVelocity2d(facexData(mph_iVelFVar, :, :, :), &
                           faceyData(mph_iVelFVar, :, :, :), &
                           facexData(mph_iRhoFVar, :, :, :), &
                           faceyData(mph_iRhoFVar, :, :, :), &
                           solnData(NRMX_VAR, :, :, :), &
                           solnData(NRMY_VAR, :, :, :), &
                           -solnData(MFLX_VAR, :, :, :), &
                           stnLimits(LOW, IAXIS), stnLimits(HIGH, IAXIS), &
                           stnLimits(LOW, JAXIS), stnLimits(HIGH, JAXIS))
#endif

#else

#ifdef MULTIPHASE_EVAPORATION
   call mph_evapVelocity3d(facexData(mph_iVelFVar, :, :, :), &
                           faceyData(mph_iVelFVar, :, :, :), &
                           facezData(mph_iVelFVar, :, :, :), &
                           facexData(mph_iRhoFVar, :, :, :), &
                           faceyData(mph_iRhoFVar, :, :, :), &
                           facezData(mph_iRhoFVar, :, :, :), &
                           solnData(NRMX_VAR, :, :, :), &
                           solnData(NRMY_VAR, :, :, :), &
                           solnData(NRMZ_VAR, :, :, :), &
                           solnData(MFLX_VAR, :, :, :), &
                           stnLimits(LOW, IAXIS), stnLimits(HIGH, IAXIS), &
                           stnLimits(LOW, JAXIS), stnLimits(HIGH, JAXIS), &
                           stnLimits(LOW, KAXIS), stnLimits(HIGH, KAXIS))
#endif

   call Stencils_advectWeno3d(solnData(HDN0_VAR, :, :, :), &
                              solnData(DFUN_VAR, :, :, :), &
                              facexData(mph_iVelFVar, :, :, :), &
                              faceyData(mph_iVelFVar, :, :, :), &
                              facezData(mph_iVelFVar, :, :, :), &
                              del(DIR_X), &
                              del(DIR_Y), &
                              del(DIR_Z), &
                              stnLimits(LOW, IAXIS), stnLimits(HIGH, IAXIS), &
                              stnLimits(LOW, JAXIS), stnLimits(HIGH, JAXIS), &
                              stnLimits(LOW, KAXIS), stnLimits(HIGH, KAXIS), &
                              center=1, facex=0, facey=0, facez=0)

#ifdef MULTIPHASE_EVAPORATION
   call mph_evapVelocity3d(facexData(mph_iVelFVar, :, :, :), &
                           faceyData(mph_iVelFVar, :, :, :), &
                           facezData(mph_iVelFVar, :, :, :), &
                           facexData(mph_iRhoFVar, :, :, :), &
                           faceyData(mph_iRhoFVar, :, :, :), &
                           facezData(mph_iRhoFVar, :, :, :), &
                           solnData(NRMX_VAR, :, :, :), &
                           solnData(NRMY_VAR, :, :, :), &
                           solnData(NRMZ_VAR, :, :, :), &
                           -solnData(MFLX_VAR, :, :, :), &
                           stnLimits(LOW, IAXIS), stnLimits(HIGH, IAXIS), &
                           stnLimits(LOW, JAXIS), stnLimits(HIGH, JAXIS), &
                           stnLimits(LOW, KAXIS), stnLimits(HIGH, KAXIS))
#endif

#endif

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)
#if NDIM == MDIM
   call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif
   call Timers_stop("Multiphase_advection")

   return
end subroutine Multiphase_advection
