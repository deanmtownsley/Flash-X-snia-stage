!!****if* source/Simulation/SimulationMain/unitTest/MoL/IMEX/Simulation_init
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
   use Simulation_data, only: sim_alpha, sim_beta, sim_A, sim_mu, sim_sigma, U_RHS

   use RuntimeParameters_interface, only: RuntimeParameters_get
   use MoL_interface, only: MoL_registerVariable

#include "Simulation.h"
#include "constants.h"

   implicit none

   call RuntimeParameters_get("sim_alpha", sim_alpha)
   call RuntimeParameters_get("sim_beta", sim_beta)
   call RuntimeParameters_get("sim_A", sim_A)
   call RuntimeParameters_get("sim_mu", sim_mu)
   call RuntimeParameters_get("sim_sigma", sim_sigma)

   call MoL_registerVariable("u", U_VAR, U_RHS)

end subroutine Simulation_init
