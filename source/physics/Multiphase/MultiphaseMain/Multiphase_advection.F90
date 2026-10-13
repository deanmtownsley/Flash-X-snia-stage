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

subroutine Multiphase_advection(solnData, facexData, faceyData, facezData,del, lo, hi)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Stencils_interface, ONLY: Stencils_advectWeno
   use mph_evapInterface, ONLY: mph_evapVelocity2d, mph_evapVelocity3d

!-----------------------------------------------------------------------------------------
   implicit none
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   real,dimension(MDIM), intent(IN) :: del
   integer, dimension(MDIM),intent(IN) :: lo, hi
   !-----------------------------------------------------------------------------------------
   integer, dimension(MDIM+1) :: face
   integer, parameter :: center=1, facex=2,facey=3,facez=4


   call Timers_start("Multiphase_advection")

   face=0

#if NDIM < MDIM

#ifdef MULTIPHASE_EVAPORATION
   call mph_evapVelocity2d(facexData(mph_iVelFVar, :, :, :), &
                           faceyData(mph_iVelFVar, :, :, :), &
                           solnData(mph_iRhoCVar, :, :, :), &
                           solnData(NRMX_VAR, :, :, :), &
                           solnData(NRMY_VAR, :, :, :), &
                           solnData(MFLX_VAR, :, :, :), &
                           lo, hi)
#endif
   face(center)=1
   call Stencils_advectWeno(solnData(HDN0_VAR, :, :, :), &
                              solnData(DFUN_VAR, :, :, :), &
                              facexData(mph_iVelFVar, :, :, :), &
                              faceyData(mph_iVelFVar, :, :, :), &
                              facezData(mph_iVelFVar, :, :, :), &
                              del, lo, hi, face)

#ifdef MULTIPHASE_EVAPORATION
   call mph_evapVelocity2d(facexData(mph_iVelFVar, :, :, :), &
                           faceyData(mph_iVelFVar, :, :, :), &
                           solnData(mph_iRhoCVar, :, :, :), &
                           solnData(NRMX_VAR, :, :, :), &
                           solnData(NRMY_VAR, :, :, :), &
                           -solnData(MFLX_VAR, :, :, :), &
                           lo, hi ) 
#endif

#else

#ifdef MULTIPHASE_EVAPORATION
   call mph_evapVelocity3d(facexData(mph_iVelFVar, :, :, :), &
                           faceyData(mph_iVelFVar, :, :, :), &
                           facezData(mph_iVelFVar, :, :, :), &
                           solnData(mph_iRhoCVar, :, :, :), &
                           solnData(NRMX_VAR, :, :, :), &
                           solnData(NRMY_VAR, :, :, :), &
                           solnData(NRMZ_VAR, :, :, :), &
                           solnData(MFLX_VAR, :, :, :), &
                           lo, hi)
#endif
   face(center)=1
   call Stencils_advectWeno(solnData(HDN0_VAR, :, :, :), &
                              solnData(DFUN_VAR, :, :, :), &
                              facexData(mph_iVelFVar, :, :, :), &
                              faceyData(mph_iVelFVar, :, :, :), &
                              facezData(mph_iVelFVar, :, :, :), &
                              del, lo, hi, face)

#ifdef MULTIPHASE_EVAPORATION
   call mph_evapVelocity3d(facexData(mph_iVelFVar, :, :, :), &
                           faceyData(mph_iVelFVar, :, :, :), &
                           facezData(mph_iVelFVar, :, :, :), &
                           solnData(mph_iRhoCVar, :, :, :), &
                           solnData(NRMX_VAR, :, :, :), &
                           solnData(NRMY_VAR, :, :, :), &
                           solnData(NRMZ_VAR, :, :, :), &
                           -solnData(MFLX_VAR, :, :, :), &
                           lo, hi)
#endif

#endif

   Call Timers_stop("Multiphase_advection")

   return
end subroutine Multiphase_advection
