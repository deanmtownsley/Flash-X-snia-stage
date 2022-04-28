!!****if* source/Simulation/SimulationMain/StreamingSineWave/Simulation_data
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
#include "Flash.h"

  !! *** Runtime Parameters *** !!
  real, save :: sim_dens_lo_i
  real, save :: sim_dens_hi_i
  real, save :: sim_temp_i

  real, save :: sim_xmin
  real, save :: sim_xmax

  !! *** Variables pertaining to this Simulation *** !!
  integer, save :: sim_meshMe, sim_globalMe
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
