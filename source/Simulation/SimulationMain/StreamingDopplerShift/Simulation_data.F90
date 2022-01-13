!!****if* source/Simulation/SimulationMain/StreamingDopplerShift/Simulation_data
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
!!  Store the simulation data for StreamingDopplerShift setup
!!
!! PARAMETERS   
!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Simulation.h"

  !! *** Variables pertaining to this Simulation *** !!

  logical, save :: sim_restart
  integer, save :: sim_meshMe, sim_globalMe
  real, save :: sim_dens_i
  real, save :: sim_temp_i
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

  character(len=MAX_STRING_LENGTH) :: sim_str_geometry
  integer, save :: sim_geometry

  !! *** Variables pertaining to this Simulation *** !!

  character(len=20), save :: sim_rad_spectrum
  character(len=1),  save :: sim_rad_direction
  real, save :: sim_velx
  real, save :: sim_s ! sqrt( (1+v/c)/(1-v/c) )

  real, parameter :: X_0 = 2.0d2
  real, parameter :: X_1 = 3.5d2
  real, parameter :: X_2 = 6.5d2
  real, parameter :: X_3 = 8.0d2
  real, parameter :: L_X = 6.0d2

end module Simulation_data
