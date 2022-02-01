!!****if* source/Simulation/SimulationMain/unitTest/Grid/Amrex/TestRefine/Simulation_data
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
!!  Store the simulation data for the Sod problem
!!
!!***

module Simulation_data
    implicit none

    integer, parameter :: MIN_REFINE_LEVEL = 1
    integer, parameter :: MAX_REFINE_LEVEL = 4
    integer, parameter :: NUM_LEVELS = MAX_REFINE_LEVEL - MIN_REFINE_LEVEL + 1 

    type blocks_t
      integer, allocatable :: blocks(:, :) 
    end type blocks_t
    
    ! Store the leaf blocks at each level
    type(blocks_t), save :: leaves(MIN_REFINE_LEVEL:MAX_REFINE_LEVEL)
end module Simulation_data

