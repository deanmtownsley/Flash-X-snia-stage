!!****if* source/Simulation/SimulationMain/CCSN_WL/Simulation_data
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
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: CCSN
!!  
!! PARAMETERS
!!
!!
!! NOTES
!!  
!!  This problem is described in, e.g.,
!!  Couch, S.M. 2013, ApJ, 765, 29
!!  Couch, S.M. 2013, ApJ, 775, 35
!!  Couch, S.M. & O'Connor, E.P. 2013, arXiv:1310.5728
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Simulation.h"
 
  !! *** Runtime Parameters *** !!

  logical, save :: sim_restart
  integer, save :: nsub
  integer, save :: sim_meshMe, sim_meshComm
  real, save    :: sim_small, sim_smallx
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  
  !! *** Variables pertaining to this Simulation *** !!
  integer, save      :: nvar_stored
  integer, parameter :: n1d_max = 10000 ! Max number of lines a file can have
  integer, save      :: n1d_total ! Actual number of lines, calculated after input
  character(len=MAX_STRING_LENGTH), save :: model_file
  real, save         :: xzn(n1d_max)
  real, save         :: model_1d(n1d_max,NUNK_VARS)
  real, save         :: sim_velMult
  character (len=4), save :: unklabels(UNK_VARS_BEGIN:UNK_VARS_END)

end module Simulation_data
