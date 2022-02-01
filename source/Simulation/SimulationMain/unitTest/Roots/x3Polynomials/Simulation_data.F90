!!****if* source/Simulation/SimulationMain/unitTest/Roots/x3Polynomials/Simulation_data
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
!!  Stores the local data for the Roots/x3Polynomial unit test.
!!  
!!***

Module Simulation_data

  implicit none

#include "constants.h"
#include "Simulation.h"

  character (len = MAX_STRING_LENGTH), save :: sim_baseName
  character (len = MAX_STRING_LENGTH), save :: sim_printName

  logical, save :: sim_printInfo

  integer, save :: sim_printUnit

  integer, parameter :: sim_numberOfx3Polynomials = 9

  real,    save :: sim_infinity

  real,    save :: sim_maxRelativeAccuracy (1:sim_numberOfx3Polynomials)
  real,    save :: sim_x3Polynomialx2Coeff (1:sim_numberOfx3Polynomials)
  real,    save :: sim_x3Polynomialx1Coeff (1:sim_numberOfx3Polynomials)
  real,    save :: sim_x3Polynomialx0Coeff (1:sim_numberOfx3Polynomials)
  real,    save :: sim_x3RelativeAccuracy  (1:sim_numberOfx3Polynomials)

  real,    save :: sim_rootsAnalytical       (1:3 , 1:2 , 1:sim_numberOfx3Polynomials)
  real,    save :: sim_rootsRelativeAccuracy (1:3 , 1:2 , 1:sim_numberOfx3Polynomials)

end module Simulation_data
