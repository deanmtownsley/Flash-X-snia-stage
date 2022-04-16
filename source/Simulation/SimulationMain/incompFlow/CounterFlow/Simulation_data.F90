!!****if* source/Simulation/SimulationMain/incompFlow/CounterFlow/Simulation_data
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
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Pool Boiling
!!
!!***

module Simulation_data

   implicit none

#include "constants.h"
#include "Simulation.h"

   real, save :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
   integer, save :: sim_meshMe

   real, save :: sim_gravX
   real, save :: sim_gravY
   real, save :: sim_gravZ

   real, save :: sim_channelDepth
   real, save :: sim_nozzleFreq
   real, save :: sim_nozzleAmp
   real, save :: sim_liqFlowRate, sim_gasFlowRate

end module Simulation_data
