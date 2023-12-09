!!****if* source/Simulation/SimulationMain/DiffusionSineWave/Simulation_init
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
!!
!! DESCRIPTION
!!
!!   Initialize general solution data for diffusion sine wave radiation
!!   problem
!!
!! ARGUMENTS
!!
!!
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Driver_interface, ONLY : Driver_abort, Driver_getComm, Driver_getMype
  use Eos_interface, ONLY : Eos_getAbarZbar
  use Logfile_interface, ONLY : Logfile_stamp
  use ProgramHeaderModule, ONLY : nE, nDOF
  use RadiationFieldsModule, ONLY : nCR, nSpecies
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"
#include "Multispecies.h"

  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction

  call Driver_getComm(MESH_COMM, sim_meshComm)
  call Driver_getMype(MESH_COMM, sim_meshMe)
  call Driver_getComm(GLOBAL_COMM, sim_globalComm)
  call Driver_getMype(GLOBAL_COMM, sim_globalMe)

  call RuntimeParameters_get('xmin', sim_xmin)
  call RuntimeParameters_get('xmax', sim_xmax)

  call RuntimeParameters_get( 'dens_i', sim_dens_i )
  call RuntimeParameters_get( 'temp_i', sim_temp_i )
  call RuntimeParameters_get( 'pres_i', sim_pres_i )
  call RuntimeParameters_get( 'vel_i' , sim_velx_i )
  
  sim_vely_i  = 0.0e0
  sim_velz_i  = 0.0e0
  sim_xn_i(:) = 0.0e0
  
#if NSPECIES > 0
  sim_xn_i(SPECIES_BEGIN) = 1.0e0
#endif

  massFraction(:) = sim_xn_i(:)

  call Eos_getAbarZbar( Ye=sim_ye_i, massFrac=massfraction )

  sim_nComp = nSpecies * nCR * nE * nDOF

  call RuntimeParameters_get( 'rt_Sigma', sim_Sigma )

  return
end subroutine Simulation_init
