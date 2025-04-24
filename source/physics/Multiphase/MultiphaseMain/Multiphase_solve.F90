!!***if* source/physics/Multiphase/MultiphaseMain/Multiphase_solve
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

subroutine Multiphase_solve(solnData, lo, hi, dt)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Stencils_interface, ONLY: Stencils_integrateEuler, Stencils_integrateAB2

   implicit none
   !----------Arugments List---------------
   real, dimension(:,:,:,:), pointer :: solnData
   real, INTENT(IN) :: dt
   integer, dimension(MDIM), intent(IN) :: lo, hi
!-----------------------------------------------------------------------------------------

   call Timers_start("Multiphase_solve")

   call Stencils_integrateEuler(solnData(DFUN_VAR, :, :, :), &
                                solnData(HDN0_VAR, :, :, :), &
                                dt, &
                                lo , hi, &
                                iSource=solnData(DFRC_VAR, :, :, :))

   call Timers_stop("Multiphase_solve")

   return
end subroutine Multiphase_solve
