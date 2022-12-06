!!****if* source/Simulation/SimulationMain/StreamingSineWave/Simulation_data
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
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for StreamingSineWave setup
!!
!! PARAMETERS   
!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Simulation.h"

  !! *** Runtime Parameters *** !!
  real, save :: sim_dens_lo_i
  real, save :: sim_dens_hi_i
  real, save :: sim_temp_i

  real, save :: sim_xmin
  real, save :: sim_xmax

  !! *** Variables pertaining to this Simulation *** !!
  integer, save :: sim_meshMe, sim_meshComm, sim_globalMe, sim_globalComm
  real, save :: sim_velx_i
  real, save :: sim_vely_i
  real, save :: sim_velz_i
  real, save :: sim_pres_i
  real, save :: sim_eint_i
  real, save :: sim_etot_i
  real, save :: sim_gamc_i
  real, save :: sim_game_i
  real, save :: sim_xn_i(SPECIES_BEGIN:SPECIES_END)
  real, save :: sim_ye_i

  real, save :: sim_nComp

end module Simulation_data
