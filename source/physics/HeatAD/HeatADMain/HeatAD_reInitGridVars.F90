!!****if* source/physics/HeatAD/HeatADMain/HeatAD_reInitGridVars
!!
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
#include "HeatAD.h"
#include "Simulation.h"

subroutine HeatAD_reInitGridVars(solnData)

   use Timers_interface, ONLY: Timers_start, Timers_stop
   use HeatAD_data

   !------------------------------------------------------------------------------------------
   implicit none
   real, pointer, dimension(:, :, :, :) :: solnData
   !------------------------------------------------------------------------------------------


   call Timers_start("HeatAD_reInitGridVars")


   
   ! DEVNOTE (10/24/2023):
   ! See accompanying changes in reInitGridVars routines
   ! in other physics units. Only including initialzation
   ! that are necessary.
#ifdef HEATAD_VARDIFFUSION
   solnData(ALPH_VAR, :,:,:) = 1.
#endif
   solnData(HTN0_VAR, :,:,:) = 0.
   solnData(TFRC_VAR, :,:,:) = 0.
   
   ! Release pointers:
   
   call Timers_stop("HeatAD_reInitGridVars")

   return
end subroutine HeatAD_reInitGridVars
