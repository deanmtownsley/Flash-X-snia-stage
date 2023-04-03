!!****if* source/Simulation/SimulationMain/unitTest/Poisson/Amrex/Simulation_init
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
!!  call Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!
!! ARGUMENTS
!!
!!   none
!!
!!***

subroutine Simulation_init()

  use Simulation_data 
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  

  call RuntimeParameters_get('xmin', sim_xMin)
  call RuntimeParameters_get('xmax', sim_xMax)
  call RuntimeParameters_get('ymin', sim_yMin)
  call RuntimeParameters_get('ymax', sim_yMax)
  call RuntimeParameters_get('zmin', sim_zMin)
  call RuntimeParameters_get('zmax', sim_zMax)
  
end subroutine Simulation_init
