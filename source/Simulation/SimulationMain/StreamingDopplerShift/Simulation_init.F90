!!****if* source/Simulation/SimulationMain/StreamingDopplerShift/Simulation_init
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
!!   Initialize general solution data for streaming doppler shift radiation
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
  use Driver_interface, ONLY : Driver_abort, Driver_getMype
  use Eos_interface, ONLY : Eos, Eos_getAbarZbar
  use Logfile_interface, ONLY : Logfile_stamp
  use ProgramHeaderModule, ONLY : nE, nDOF, swE
  use RadiationFieldsModule, ONLY : nCR, nSpecies
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstantsModule, ONLY : SpeedOfLightCGS
  use UnitsModule, ONLY : Centimeter, Second, Gram
  use KindModule, ONLY : One

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"
#include "Multispecies.h"

  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction
  real, dimension(EOS_NUM) :: eosData
  real, parameter :: lightspeed = SpeedOfLightCGS
  real, parameter :: conv_dens = Gram / Centimeter**3

  call Logfile_stamp( 'Entering simulation initialization' , '[Simulation_init]')
  call Driver_getMype(MESH_COMM, sim_meshMe)
  call Driver_getMype(GLOBAL_COMM, sim_globalMe)

  call RuntimeParameters_get( 'restart', sim_restart)
  call RuntimeParameters_get("geometry",sim_str_geometry)
  call RuntimeParameters_mapStrToInt(sim_str_geometry, sim_geometry)

  call RuntimeParameters_get("sim_rad_spectrum", sim_rad_spectrum)
  call RuntimeParameters_get("sim_rad_direction", sim_rad_direction)
  call RuntimeParameters_get("sim_velx", sim_velx)

  sim_s = SQRT( ( One + sim_velx ) / ( One - sim_velx ) )

  call RuntimeParameters_get('dens_i', sim_dens_i)
  call RuntimeParameters_get('temp_i', sim_temp_i)

  sim_xn_i(:) = 0.0e0
  sim_xn_i(SPECIES_BEGIN) = 1.0e0

  massFraction(:) = sim_xn_i(:)

  eosData(EOS_DENS) = sim_dens_i * conv_dens
  eosData(EOS_TEMP) = sim_temp_i

  call Eos(MODE_DENS_TEMP,1,eosData,massFraction)

  sim_velx_i = 0.0e0
  sim_vely_i = 0.0e0
  sim_velz_i = 0.0e0

  sim_pres_i = eosData(EOS_PRES)
  sim_eint_i = eosData(EOS_EINT)
  sim_gamc_i = eosData(EOS_GAMC)
  sim_game_i = sim_pres_i/(sim_eint_i*sim_dens_i) + 1.0e0
  call Eos_getAbarZbar(Ye=sim_ye_i,massFrac=massfraction)

  if( sim_rad_direction == "X" ) sim_velx_i = sim_velx * lightspeed
  if( sim_rad_direction == "Y" ) sim_vely_i = sim_velx * lightspeed
  if( sim_rad_direction == "Z" ) sim_velz_i = sim_velx * lightspeed

  sim_etot_i = sim_eint_i + 0.5*(sim_velx_i**2 + sim_vely_i**2 + sim_velz_i**2)

  return

end subroutine Simulation_init
