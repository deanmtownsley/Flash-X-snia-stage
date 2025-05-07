!!***if* source/physics/Multiphase/MultiphaseEvap/Multiphase_thermalForcing
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

#include "constants.h"
#include "Multiphase.h"
#include "Simulation.h"

subroutine Multiphase_thermalForcing(solnData, del, logc, higc)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use mph_evapInterface, ONLY: mph_tempGfm2d, mph_tempGfm3d

!------------------------------------------------------------------------------------------------
   implicit none

   integer, dimension(MDIM),intent(IN) :: logc, higc 
   real, pointer, dimension(:, :, :, :) :: solnData
   real,dimension(MDIM), intent(IN) :: del

   call Timers_start("Multiphase_thermalForcing")

#if NDIM < MDIM

   call mph_tempGfm2d(solnData(DFUN_VAR, :, :, :), &
                      solnData(NRMX_VAR, :, :, :), &
                      solnData(NRMY_VAR, :, :, :), &
                      (mph_invReynolds/mph_Prandtl)*solnData(mph_iAlphaCVar, :, :, :), &
                      solnData(mph_iTempVar, :, :, :), &
                      solnData(mph_iTempFrcVar, :, :, :), &
                      solnData(HFLQ_VAR, :, :, :), &
                      solnData(HFGS_VAR, :, :, :), &
                      mph_Tsat, &
                      del(DIR_X), del(DIR_Y), &
                      logc(IAXIS),higc(IAXIS), &
                      logc(JAXIS),higc(JAXIS), tol=mph_tempTol)

#else

   call mph_tempGfm3d(solnData(DFUN_VAR, :, :, :), &
                      solnData(NRMX_VAR, :, :, :), &
                      solnData(NRMY_VAR, :, :, :), &
                      solnData(NRMZ_VAR, :, :, :), &
                      (mph_invReynolds/mph_Prandtl)*solnData(mph_iAlphaCVar, :, :, :), &
                      solnData(mph_iTempVar, :, :, :), &
                      solnData(mph_iTempFrcVar, :, :, :), &
                      solnData(HFLQ_VAR, :, :, :), &
                      solnData(HFGS_VAR, :, :, :), &
                      mph_Tsat, &
                      del(DIR_X), del(DIR_Y), del(DIR_Z), &
                      logc(IAXIS),higc(IAXIS), &
                      logc(JAXIS),higc(JAXIS), &
                      logc(KAXIS),higc(KAXIS), tol=mph_tempTol)

#endif


   call Timers_stop("Multiphase_thermalForcing")

   return

end subroutine Multiphase_thermalForcing
