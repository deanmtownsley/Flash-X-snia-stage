!!***if* source/physics/HeatAD/HeatADMain/varDiffusion/HeatAD_diffusion
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

subroutine HeatAD_diffusion(solnData, del, lo, hi)

   use HeatAD_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Stencils_interface, ONLY: Stencils_diffusion

!--------------------------------------------------------------------------------------------
   implicit none
   real, pointer, dimension(:, :, :, :) :: solnData
   real,dimension(MDIM),intent(IN) ::  del
   integer, dimension(MDIM),intent(IN) :: lo, hi
   real :: diffusion_coeff

!---------------------------------------------------------------------------------------------

   call Timers_start("HeatAD_diffusion")

   diffusion_coeff = ht_invReynolds/ht_Prandtl

  
   call Stencils_diffusion(solnData(HTN0_VAR, :, :, :), &
                             solnData(TEMP_VAR, :, :, :), &
                             del,&
                             diffusion_coeff*solnData(ALPH_VAR, :, :, :), &
                             lo, hi)


   call Timers_stop("HeatAD_diffusion")

end subroutine HeatAD_diffusion
