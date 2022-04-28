!!****if* source/Simulation/SimulationMain/unitTest/Eos/Simulation_init
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
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  call Simulation_init( )
!!
!!
!! DESCRIPTION
!!  Initializes all the parameters needed the Eos unit test
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data, ONLY : sim_xmin,sim_xmax,sim_ymin,sim_ymax,&
                              sim_zmin,sim_zmax,sim_smallx,sim_smallE,&
                              sim_densMin,sim_tempMin,sim_yeMin,sim_presMin,&
                              sim_densMax, sim_tempMax, sim_yeMax, sim_presMax, &
                              sim_initialMass
  use Simulation_data, ONLY : sim_meshMe, sim_debug
  use Driver_interface, ONLY : Driver_abort, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

# include "Simulation.h"
# include "constants.h"

  

  integer :: lrefine_max, factor, nblockx,nblocky,nblockz

  call Driver_getMype(MESH_COMM,   sim_meshMe)
  call RuntimeParameters_get('sim_debug',sim_debug)

  call eos_initTest
  
  call RuntimeParameters_get( 'xmin', sim_xmin)
  call RuntimeParameters_get( 'xmax', sim_xmax)
  
  call RuntimeParameters_get( 'ymin', sim_ymin)
  call RuntimeParameters_get( 'ymax', sim_ymax)
  
  call RuntimeParameters_get( 'zmin', sim_zmin)
  call RuntimeParameters_get( 'zmax', sim_zmax)
  
  call RuntimeParameters_get( 'smallx', sim_smallx)
  call RuntimeParameters_get( 'smallE', sim_smallE)

! sim_initialMass must be less than NSPECIES
  call RuntimeParameters_get( 'sim_initialMass', sim_initialMass)
  if (sim_initialMass .GT. NSPECIES)                        &
       call Driver_abort('[Simulation_init] sim_initialMass must be less than NSPECIES')


  call RuntimeParameters_get( 'sim_densMin', sim_densMin)
  call RuntimeParameters_get( 'sim_densMax', sim_densMax)
  
  call RuntimeParameters_get( 'sim_tempMin', sim_tempMin)
  call RuntimeParameters_get( 'sim_tempMax', sim_tempMax)

  call RuntimeParameters_get( 'sim_presMin', sim_presMin)
  call RuntimeParameters_get( 'sim_presMax', sim_presMax)
  
  call RuntimeParameters_get( 'sim_yeMin', sim_yeMin)
  call RuntimeParameters_get( 'sim_yeMax', sim_yeMax)

end subroutine Simulation_init
