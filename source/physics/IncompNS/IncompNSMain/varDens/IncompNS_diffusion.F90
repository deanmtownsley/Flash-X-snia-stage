!!****if* source/physics/IncompNS/IncompNSMain/varDens/IncompNS_diffusion
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
!!***
!!REORDER(4): face[xyz]Data
!!REORDER(4): solnData

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

subroutine IncompNS_diffusion(solnData, facexData, faceyData, facezData, del, lo, hi)

   use ins_interface, ONLY: ins_diffusion2d_vardens, ins_diffusion3d_vardens
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use IncompNS_data

!------------------------------------------------------------------------------------------
   implicit none
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   real, dimension(MDIM), intent(in) :: del
   integer,dimension(MDIM), intent(in) :: lo,hi

!------------------------------------------------------------------------------------------
   !
   call Timers_start("IncompNS_diffusion")
   !

#if NDIM == 3
   ! compute RHS of momentum equation
   call ins_diffusion3d_vardens(facexData(VELC_FACE_VAR, :, :, :), &
                                faceyData(VELC_FACE_VAR, :, :, :), &
                                facezData(VELC_FACE_VAR, :, :, :), &
                                solnData(TVIS_VAR, :, :, :), &
                                ins_invReynolds, lo, hi, del, &
                                facexData(HVN0_FACE_VAR, :, :, :), &
                                faceyData(HVN0_FACE_VAR, :, :, :), &
                                facezData(HVN0_FACE_VAR, :, :, :), &
                                solnData(VISC_VAR, :, :, :), &
                                facexData(RHOF_FACE_VAR, :, :, :), &
                                faceyData(RHOF_FACE_VAR, :, :, :), &
                                facezData(RHOF_FACE_VAR, :, :, :))

#elif NDIM ==2
   ! compute RHS of momentum equation
   call ins_diffusion2d_vardens(facexData(VELC_FACE_VAR, :, :, :), &
                                faceyData(VELC_FACE_VAR, :, :, :), &
                                ins_invReynolds,lo, hi, del, &
                                facexData(HVN0_FACE_VAR, :, :, :), &
                                faceyData(HVN0_FACE_VAR, :, :, :), &
                                solnData(VISC_VAR, :, :, :), &
                                facexData(RHOF_FACE_VAR, :, :, :), &
                                faceyData(RHOF_FACE_VAR, :, :, :))

#endif

   call Timers_stop("IncompNS_diffusion")

   return
end subroutine IncompNS_diffusion
