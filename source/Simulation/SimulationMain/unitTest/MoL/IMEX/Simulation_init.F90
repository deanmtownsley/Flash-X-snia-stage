!!****if* source/Simulation/SimulationMain/unitTest/MoL/IMEX/Simulation_init
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  call Simulation_init()
!!
!!
!! DESCRIPTION
!!  Initializes all the parameters needed for a particular simulation
!!
!! ARGUMENTS
!!
!!
!!***
subroutine Simulation_init()
   use Simulation_data, only: sim_alpha, sim_beta, sim_epsilon, &
                              sim_lambdaF, sim_lambdaS, sim_k, &
                              U_RHS, V_RHS

   use RuntimeParameters_interface, only: RuntimeParameters_get
   use MoL_interface, only: MoL_registerVariable, MoL_getRHSIndex

#include "Simulation.h"
#include "constants.h"

   implicit none

   call RuntimeParameters_get("sim_alpha", sim_alpha)
   call RuntimeParameters_get("sim_beta", sim_beta)
   call RuntimeParameters_get("sim_epsilon", sim_epsilon)
   call RuntimeParameters_get("sim_lambdaF", sim_lambdaF)
   call RuntimeParameters_get("sim_lambdaS", sim_lambdaS)

   call RuntimeParameters_get("sim_k", sim_k)

   call MoL_registerVariable("u", U_VAR, U_RHS)
   call MoL_registerVariable("v", V_VAR, V_RHS)

end subroutine Simulation_init
