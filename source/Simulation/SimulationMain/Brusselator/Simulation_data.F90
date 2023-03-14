!!****if* source/Simulation/SimulationMain/Brusselator/Simulation_data
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
!!  NAME
!!    Simulation_data
!!
!!  SYNOPSIS
!!    use Simulation_data
!!
!!  DESCRIPTION
!!
!!    Stores the local data for Simulation setup
!!
!!***
module Simulation_data

   implicit none

   real, save :: sim_a, sim_b
   real, save :: sim_epsilon, sim_alpha, sim_rho

   integer, save :: sim_k

   integer, save :: U_RHS, V_RHS, W_RHS

end module Simulation_data
