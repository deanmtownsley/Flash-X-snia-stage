!!***if* source/physics/HeatAD/HeatADMain/HeatAD_solve
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

subroutine HeatAD_solve(solndata, del, lo, hi, dt)

   use HeatAD_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY:  Driver_abort
   use Stencils_interface, ONLY: Stencils_integrateEuler, Stencils_integrateAB2

   implicit none

!--------------------------------------------------------------------------------------------
   real, pointer, dimension(:, :, :, :) :: solnData
   real, INTENT(IN) :: dt
   real,dimension(MDIM),intent(IN) ::  del
   integer, dimension(MDIM), intent(IN) :: lo, hi
   real :: diffusion_coeff

!---------------------------------------------------------------------------------------------
   call Timers_start("HeatAD_solve")


   if (ht_intSchm /= 1 .and. ht_intSchm /= 2) then
      call Driver_abort("[HeatAD_solve] ht_intSchm should be 1 or 2")
   end if

   diffusion_coeff = ht_invReynolds/ht_Prandtl


   if (ht_intSchm == 1) then
      call Stencils_integrateEuler(solnData(TEMP_VAR, :, :, :), &
                                   solnData(HTN0_VAR, :, :, :), &
                                   dt, lo, hi, &
                                   iSource=solnData(TFRC_VAR, :, :, :))

   else if (ht_intSchm == 2) then
      call Stencils_integrateAB2(solnData(TEMP_VAR, :, :, :), &
                                 solnData(HTN0_VAR, :, :, :), &
                                 solnData(HTN1_VAR, :, :, :), &
                                 dt, lo, hi, &
                                 iSource=solnData(TFRC_VAR, :, :, :))

      solnData(HTN1_VAR, :, :, :) = solnData(HTN0_VAR, :, :, :)

   end if


   call Timers_stop("HeatAD_solve")

end subroutine HeatAD_solve
