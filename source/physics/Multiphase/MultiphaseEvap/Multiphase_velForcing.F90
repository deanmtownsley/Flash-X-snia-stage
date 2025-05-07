!!***if* source/physics/Multiphase/MultiphaseEvap/Multiphase_velForcing
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
!!REORDER(4): solnData,face[xyz]Data

#include "constants.h"
#include "Multiphase.h"
#include "Simulation.h"

subroutine Multiphase_velForcing(solnData, facexData, faceyData, facezData, del, lo, hi, dt)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use mph_evapInterface, ONLY: mph_evapVelForcing2d, mph_evapVelForcing3d

   implicit none
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   real, dimension(MDIM), intent(in) :: del
   integer,dimension(MDIM), intent(in) :: lo,hi
   real, INTENT(IN) :: dt

   
!------------------------------------------------------------------------------------------------
   integer :: ierr, i, j, k
   integer TA(2), count_rate
   real*8 ET
!------------------------------------------------------------------------------------------------
   
#if NDIM < MDIM
   call mph_evapVelForcing2d(facexData(mph_iVelFVar, :, :, :), &
                             faceyData(mph_iVelFVar, :, :, :), &
                             facexData(mph_iRhoFVar, :, :, :), &
                             faceyData(mph_iRhoFVar, :, :, :), &
                             solnData(mph_iRhoCVar, :, :, :), &
                             solnData(mph_iMuCVar, :, :, :), &
                             solnData(NRMX_VAR, :, :, :), &
                             solnData(NRMY_VAR, :, :, :), &
                             solnData(MFLX_VAR, :, :, :), &
                             mph_invReynolds, dt, del, lo,hi)

#else

   call mph_evapVelForcing3d(facexData(mph_iVelFVar, :, :, :), &
                             faceyData(mph_iVelFVar, :, :, :), &
                             facezData(mph_iVelFVar, :, :, :), &
                             facexData(mph_iRhoFVar, :, :, :), &
                             faceyData(mph_iRhoFVar, :, :, :), &
                             facezData(mph_iRhoFVar, :, :, :), &
                             solnData(mph_iRhoCVar, :, :, :), &
                             solnData(mph_iMuCVar, :, :, :), &
                             solnData(NRMX_VAR, :, :, :), &
                             solnData(NRMY_VAR, :, :, :), &
                             solnData(NRMZ_VAR, :, :, :), &
                             solnData(MFLX_VAR, :, :, :), &
                             mph_invReynolds, dt, del, lo, hi)

#endif

   return

end subroutine Multiphase_velForcing
