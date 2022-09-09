!!****if* source/numericalTools/MoL/MoLMain/MR/ml_data
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
!!  NAME
!!    ml_data
!!
!!  SYNOPSIS
!!    use ml_data
!!
!!  DESCRIPTION
!!    Stores data for the multi-rate implicit-explicit integrator
!!
!!  NOTES
!!
!!***
module ml_data

   implicit none

   character(len=:), allocatable, save :: ml_slowMethod
   character(len=:), allocatable, save :: ml_fastMethod

   integer, save :: ml_slowOrder, ml_fastOrder
   integer, save :: ml_nstages_slow, ml_nstages_fast
   integer, save :: ml_nsubcycle

   ! Note: indexing variables below do not contain an ml_ to maintain a consistent notation
   ! with the pre-processor defined variable indexing in Simulation.h and elsewhere

   ! RHS indexinng
   integer, dimension(:), allocatable, save :: FE, FI, FF

   ! Initial state index for fast evolution
   integer, save :: FAST_INITIAL

   ! Butcher tableau for slow evolution
   integer, save :: ml_kmax
   real, dimension(:, :, :), allocatable, save :: ml_gamK, ml_wK
   real, dimension(:, :), allocatable, save :: ml_gamBar, ml_wBar
   real, dimension(:), allocatable, save :: ml_cS

   ! Butcher tableau for fast evolution
   real, dimension(:, :), allocatable, save :: ml_AF
   real, dimension(:), allocatable, save :: ml_bF, ml_cF

end module ml_data
