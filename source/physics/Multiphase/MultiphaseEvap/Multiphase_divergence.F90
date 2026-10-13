!!***if* source/physics/Multiphase/MultiphaseEvap/Multiphase_divergence
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

subroutine Multiphase_divergence(solnData, facexData, faceyData, facezData, del, lo, hi)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Stencils_interface, ONLY: Stencils_cnt_advectUpwind2d, Stencils_cnt_advectUpwind3d
   use mph_evapInterface, ONLY: mph_evapDivergence2d, mph_evapDivergence3d

!------------------------------------------------------------------------------------------------
   implicit none
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   real,dimension(MDIM), intent(IN) :: del
   integer, dimension(MDIM),intent(IN) :: lo, hi
 !------------------------------------------------------------------------------------------------

   call Timers_start("Multiphase_divergence")

#if NDIM < MDIM
   call mph_evapDivergence2d(solnData(mph_iDivCvar, :, :, :), &
                             solnData(mph_iRhoCVar, :, :, :), &
                             solnData(NRMX_VAR, :, :, :), &
                             solnData(NRMY_VAR, :, :, :), &
                             solnData(MFLX_VAR, :, :, :), &
                             del(DIR_X), del(DIR_Y), &
                             lo(IAXIS), hi(IAXIS), &
                             lo(JAXIS), hi(JAXIS))

#else

   call mph_evapDivergence3d(solnData(mph_iDivCvar, :, :, :), &
                             solnData(mph_iRhoCVar, :, :, :), &
                             solnData(NRMX_VAR, :, :, :), &
                             solnData(NRMY_VAR, :, :, :), &
                             solnData(NRMZ_VAR, :, :, :), &
                             solnData(MFLX_VAR, :, :, :), &
                             del(DIR_X), del(DIR_Y), del(DIR_Z), &
                             lo(IAXIS), hi(IAXIS), &
                             lo(JAXIS), hi(JAXIS), &
                             lo(KAXIS), hi(KAXIS))

#endif
   call Timers_stop("Multiphase_divergence")
   return

end subroutine Multiphase_divergence
