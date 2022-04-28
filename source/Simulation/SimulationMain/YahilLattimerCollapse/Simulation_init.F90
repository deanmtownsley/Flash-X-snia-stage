!!****if* source/Simulation/SimulationMain/YahilLattimerCollapse/Simulation_init
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
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Yahil-Lattimer Collapse
!!  problem.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!  sim_profFileName   Name of file from which to read a 1D 
!!                     YahilLattimerCollapse profile for the initial condition
!!
!!***

subroutine Simulation_init()

  use Simulation_data 
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface,           ONLY : Logfile_stamp
  use ut_generalInterface,         ONLY : ut_getFreeFileUnit

  implicit none

#include "Simulation.h"
#include "constants.h"

  logical :: threadBlockListBuild, threadWithinBlockBuild  
  real    :: vctr

  call Driver_getMype(MESH_COMM,   sim_meshMe)
  call Driver_getMype(GLOBAL_COMM, sim_globalMe)
  call Driver_getNumProcs(GLOBAL_COMM, sim_globalNumProcs)

  call RuntimeParameters_get('sim_pCentral', sim_pCentral)
  call RuntimeParameters_get('sim_rhoCentral', sim_rhoCentral)
  call RuntimeParameters_get('sim_gammaInitial', sim_gammaInitial)
  call RuntimeParameters_get('sim_collapsetime', sim_collapsetime)
  call RuntimeParameters_get('sim_maxDens', sim_maxDens)
  call RuntimeParameters_get('sim_xctr',sim_xCenter)
  call RuntimeParameters_get('sim_yctr',sim_yCenter)
  call RuntimeParameters_get('tinitial',sim_tInitial)
  call RuntimeParameters_get('sim_zctr',sim_zCenter)
  call RuntimeParameters_get('nsub',sim_nSubZones)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smallt', sim_smallT)
  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('ymin',sim_yMin)
  call RuntimeParameters_get('zmin',sim_zMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  call RuntimeParameters_get('ymax',sim_yMax)
  call RuntimeParameters_get('zmax',sim_zMax)
  call RuntimeParameters_get('sim_profFileName',sim_profFileName)

  sim_useProfileFromFile = .FALSE.
  if (sim_nProfile > 1) then
     sim_useProfileFromFile = (trim(sim_profFileName) .NE. "/dev/null")
  end if

  sim_inSubZones = 1./real(sim_nSubZones)
  sim_inSubzm1   = 1./real(max(sim_nSubZones-1,1))
  sim_inszd      = sim_inSubZones**NDIM

  if (sim_useProfileFromFile) then
     call sim_readProfile()
  end if

  sim_postBounce = .false.  

end subroutine Simulation_init
