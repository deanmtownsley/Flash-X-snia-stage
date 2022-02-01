!!****if* source/Simulation/SimulationMain/HydroStatic/Simulation_data
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
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Hydrostatic
!!  
!!
!!
!!***
module Simulation_data
  implicit none
#include "constants.h"
  real, save :: sim_gamma, sim_smallX
  integer, save :: sim_gravVector(MDIM)
  real, save :: sim_GravConst
  integer, save :: sim_gravDirec

  real, save :: sim_xyzRef, sim_presRef, sim_densRef, sim_tempRef
  real, save :: sim_molarMass, sim_gasconstant
integer, save :: sim_meshMe
end module Simulation_data
