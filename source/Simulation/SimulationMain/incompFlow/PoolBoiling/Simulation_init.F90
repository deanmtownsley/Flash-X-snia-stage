!!****if* source/Simulation/SimulationMain/incompFlow/PoolBoiling/Simulation_init
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
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! ARGUMENTS
!!
!!   none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for INS-isotropic turbulence problem.
!!
!!***

#include "constants.h"

subroutine Simulation_init()

  use Simulation_data
  use Driver_interface,            ONLY : Driver_getMype,Driver_abort
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use sim_heaterInterface,         ONLY : sim_heaterInit
  use sim_outflowInterface,        ONLY : sim_outflowInit

  implicit none

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  call RuntimeParameters_get('zmin',    sim_zmin)
  call RuntimeParameters_get('zmax',    sim_zmax)

  call sim_heaterInit()
  call sim_outflowInit()

end subroutine Simulation_init
