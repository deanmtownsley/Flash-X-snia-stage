!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson3/Simulation_data
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
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup:  MacLaurin
!!  
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"
  
  !! *** Runtime Parameters *** !!
  real, save    :: sim_passTolerance !! percentage error allowed in pass/fail
  
  real, save    :: sim_e  !! eccentricity
  real, save    :: sim_gamma, sim_density, sim_Omega1, sim_a1
  real, save    :: sim_xctr, sim_yctr, sim_zctr
  real, save    :: sim_smallX, sim_smallE, sim_smallRho, sim_smallP
  real, save    :: sim_Newton, sim_pi

  integer, save :: sim_nsubzones, sim_initGeometry
!!  integer, parameter :: sim_geom2DAxisymmetric=1, sim_geom3DCartesian=2
!!  should use constants.h geometry = CYLINDRICAL,  CARTESIAN

  !! *** Auxiliary Variables - Introduced in Simulation_init *** !!

  real, save    :: sim_nsubinv, sim_a3, sim_a1inv, sim_a3inv
  real, save    :: sim_Omega2, sim_Pconst


integer, save :: sim_meshMe
end module Simulation_data
