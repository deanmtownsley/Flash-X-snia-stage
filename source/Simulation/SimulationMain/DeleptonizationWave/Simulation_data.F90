!!****if* source/Simulation/SimulationMain/DeleptonizationWave/Simulation_data
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
!!  Store the simulation data for DeleptonizationWave setup
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
  integer, save :: sim_meshMe, sim_meshComm, sim_globalMe, sim_globalComm

  real, save :: sim_nComp

  character(len=MAX_STRING_LENGTH) :: sim_str_geometry
  integer, save :: sim_geometry

 
  !! *** Variables pertaining to this Simulation *** !!
  ! --- Density Profile ---
  real, parameter :: MinD = 1.0d08 !* Gram / Centimeter**3
  real, parameter :: MaxD = 4.0d14 !* Gram / Centimeter**3
  real, parameter :: R_D  = 2.0d06 !* centimeter
  real, parameter :: H_D  = 1.0d06 !* centimeter
  ! --- Temperature Profile ---
  real, parameter :: MinT = 5.0d09 !* Kelvin
  real, parameter :: MaxT = 2.6d11 !* Kelvin
  real, parameter :: R_T  = 2.5d06 !* centimeter
  real, parameter :: H_T  = 2.0d06 !* centimeter
  ! --- Electron Fraction Profile ---
  real, parameter :: MinY = 3.0d-1
  real, parameter :: MaxY = 4.6d-1
  real, parameter :: R_Y  = 4.5d06 !* centimeter
  real, parameter :: H_Y  = 1.0d06 !* centimeter  

  ! uniform initial conditions
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
  real, save :: sim_mu_i

  ! profile use
  logical, save :: sim_use_model
  character(len=MAX_STRING_LENGTH), save :: sim_model_file
  integer, save :: sim_rad_option ! 0:FD, 1:Analytic, 2:chimera profile
                                  ! 3:boltztran profile
  real, save :: sim_rintSwitch

  integer, save :: nvar_stored
  integer, save :: n1d_max
  integer, save :: n1d_nE
  real, allocatable, save :: xzn(:)  ! fluid radius arrary
  real, allocatable, save :: xznrad(:) ! rad radius array
  real, allocatable, save :: ezn(:)  ! energy arrary
  real, allocatable, save :: model_1d(:,:)
  real, allocatable, save :: model_1d_rad(:,:,:,:) ! e,r,icr,is
  real, allocatable, save :: D_Nu_P(:,:,:) ! e,r,is
  real, allocatable, save :: I1_Nu_P(:,:,:) ! e,r,is

end module Simulation_data
