!!****if* source/Simulation/SimulationMain/unitTest/Eos/Simulation_data
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
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!   
!!
!!***

module Simulation_data
  implicit none
  real, save :: sim_xmin, sim_xmax, sim_ymin, sim_ymax, sim_zmin, sim_zmax

  real, save :: sim_smallx, sim_smallE

  integer, save ::  sim_initialMass

  real, save :: sim_densMin, sim_tempMin, sim_yeMin, sim_presMin
  real, save :: sim_densMax, sim_tempMax, sim_yeMax, sim_presMax

  integer, save :: sim_meshMe
  logical, save :: sim_debug
end module Simulation_data
