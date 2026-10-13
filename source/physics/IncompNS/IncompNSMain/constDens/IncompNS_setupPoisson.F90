!!****if* source/physics/IncompNS/IncompNSMain/constdens/IncompNS_setupPoisson
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
!!REORDER(4): face[xyz]Data
!!REORDER(4): solnData

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

subroutine IncompNS_setupPoisson(solnData, facexData, faceyData, facezData, del,lo, hi, dt)

   use ins_interface, ONLY: ins_setupPoissonRhs_constdens
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use IncompNS_data

   implicit none
   real, INTENT(IN) :: dt
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   real,dimension(MDIM),intent(IN) :: del
   integer, dimension(MDIM) :: lo, hi


   call Timers_start("IncompNS_setupPoisson")

   ! Poisson RHS source vector
   call ins_setupPoissonRhs_constdens(solnData(DUST_VAR, :, :, :), dt)

   ! Release pointers:

   call Timers_stop("IncompNS_setupPoisson")

   return
end subroutine IncompNS_setupPoisson
