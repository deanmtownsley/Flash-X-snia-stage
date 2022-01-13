!!****if* source/Simulation/SimulationMain/incompFlow/PoolBoiling/localAPI/outflow/sim_outflowData
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!! NAME
!!
!!  sim_outflowData
!!
!! SYNOPSIS
!!
!!  use sim_outflowData
!!
!!***

#include "constants.h"
#include "Simulation.h"

module sim_outflowData

    implicit none

    real, save :: sim_outflowVel(LOW:HIGH,MDIM)
    real, save :: sim_gravX
    real, save :: sim_gravY
    real, save :: sim_gravZ
    real, save :: sim_lsSink
    real, save :: sim_outflowRegion

end module sim_outflowData
