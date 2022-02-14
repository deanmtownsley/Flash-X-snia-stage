!!****if* source/Simulation/SimulationMain/DustCollapse/Simulation_data
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
!!  Store the simulation data for the DustCollapse problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!   sim_initRad                  Initial radius of cloud
!!   sim_initDens                 Initial density of cloud
!!   sim_tAmbient                 Initial ambient temperature (everywhere)
!!   sim_iCtr,sim_jCtr, sim_kCtr  Coordinates of the center of the cloud
!!
!!***
Module Simulation_data
  integer,parameter :: N_prof = 1000
  real, save    :: sim_imin, sim_imax, sim_jmin, sim_jmax, sim_kmin, sim_kmax
  real, save    :: sim_smalle, sim_smallp, sim_initRad, sim_gamma, sim_tAmbient
  real, save    :: sim_initDens, sim_smlrho, sim_ictr, sim_jctr, sim_kctr
  real, save    :: sim_presFrac
  real,save, dimension(N_prof):: sim_rProf, sim_rhoProf, sim_pProf,sim_vProf
  integer, save :: sim_meshMe
end Module Simulation_data

