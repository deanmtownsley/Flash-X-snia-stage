!!****if* source/Simulation/SimulationMain/Brusselator/Simulation_init
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
   use Simulation_data, only: sim_a, sim_b, sim_alpha, &
                              sim_epsilon, sim_rho, sim_k, &
                              U_RHS, V_RHS, W_RHS

   use RuntimeParameters_interface, only: RuntimeParameters_get
   use MoL_interface, only: MoL_registerVariable, MoL_getRHSIndex

#include "Simulation.h"
#include "constants.h"

   implicit none

   call RuntimeParameters_get("sim_a", sim_a)
   call RuntimeParameters_get("sim_b", sim_b)

   call RuntimeParameters_get("sim_alpha", sim_alpha)
   call RuntimeParameters_get("sim_epsilon", sim_epsilon)
   call RuntimeParameters_get("sim_rho", sim_rho)

   call RuntimeParameters_get("sim_k", sim_k)

   call MoL_registerVariable("u", U_VAR, U_RHS)
   call MoL_registerVariable("v", V_VAR, V_RHS)
   call MoL_registerVariable("w", W_VAR, W_RHS)

   ! Sanity check
   U_RHS = MoL_getRHSIndex(U_VAR)
   V_RHS = MoL_getRHSIndex(V_VAR)
   W_RHS = MoL_getRHSIndex(W_VAR)

end subroutine Simulation_init
