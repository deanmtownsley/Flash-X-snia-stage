!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/3Dcircle/Simulation_data
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
!!  Stores the local data for the Runge Kutta unit test.
!!  
!!***

Module Simulation_data

  implicit none

#include "constants.h"
#include "Simulation.h"

  character (len = MAX_STRING_LENGTH), save :: sim_RungeKuttaMethod

  integer, save :: sim_numberOfCircles
  integer, save :: sim_numberOfRungeKuttaSteps

  real,    save :: sim_errorFraction
  real,    save :: sim_radius
  real,    save :: sim_speed
  real,    save :: sim_stepSize

  real,    save :: sim_rx0
  real,    save :: sim_ry0
  real,    save :: sim_rz0

  real,    save :: sim_vx0
  real,    save :: sim_vy0
  real,    save :: sim_vz0

end module Simulation_data
