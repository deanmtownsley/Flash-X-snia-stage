!!****f* source/Simulation/SimulationMain/incompFlow/PoolBoiling/Simulation_initRestart
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
!!  Simulation_initRestart
!!
!! SYNOPSIS
!!  Simulation_initRestart()
!!
!! DESCRIPTION
!!  This is where the user should place code for a setup that needs to adjust
!!  data on a restart, particularly if grid data, grid metadata or particle 
!!  data needs to be changed on restarting.
!!
!! ARGUMENTS
!!
!!***

subroutine Simulation_initRestart()


  use Simulation_data
  use Driver_interface,            ONLY : Driver_getMype,Driver_abort
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use sim_heaterInterface,         ONLY : sim_heaterInit
  use sim_outflowInterface,        ONLY : sim_outflowInit

#include "constants.h"

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

  !if(sim_chkptSiteNum .le. sim_nucSiteNum) sim_nucTimeStampAll(1:sim_chkptSiteNum) = sim_chkptTimeStampAll
  !if(sim_chkptSiteNum .gt. sim_nucSiteNum) sim_nucTimeStampAll = sim_chkptTimeStampAll(1:sim_nucSiteNum)
  !
  !deallocate(sim_chkptTimeStampAll)

end subroutine Simulation_initRestart
