!!****if* source/Simulation/SimulationForcing/incompFlow/Inlet/sim_inletData
!!
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
!!  sim_inletData
!!
!! SYNOPSIS
!!
!!  use sim_inletData
!!
!!***

#include "constants.h"
#include "Simulation.h"

module sim_inletData

   implicit none

   integer, save :: sim_inletFlag(LOW:HIGH, MDIM)

   real, save :: sim_inletBuffer
   real, save :: sim_inletGrowthRate
   real, save :: sim_inletSink

end module sim_inletData
