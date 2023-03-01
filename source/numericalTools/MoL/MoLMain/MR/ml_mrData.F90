!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief Main data module for MoL's multi-rate integrator

!> @ingroup MoLMR
!! Stores data for the multi-rate (MR) integrator
module ml_mrData

   implicit none

   !> @name Multi-Rate Method Information
   !! @{

   !> Name of the "slow" method
   character(len=:), allocatable, save :: ml_slowMethod

   !> Name of the "fast" method
   character(len=:), allocatable, save :: ml_fastMethod

   !> Order of the "slow" method
   integer, save :: ml_slowOrder

   !> Order of the "fast" method
   integer, save :: ml_fastOrder

   !> Number of "slow" stages
   integer, save :: ml_nstages_slow

   !> Number of "fast" stages
   integer, save :: ml_nstages_fast

   !> Number of steps taken during each fast stage
   integer, save :: ml_nsubcycle
   !> @}

   !> @name Indexing
   !! @{
   !! @note Indexing variables below do not contain an ml_ to maintain a consistent notation
   !!       with the pre-processor defined variable indexing in Simulation.h and elsewhere

   !> Explicit RHS indices
   integer, dimension(:), allocatable, save :: FE

   !> Implicit RHS indices
   integer, dimension(:), allocatable, save :: FI

   !> Fast RHS indices
   integer, dimension(:), allocatable, save :: FF

   !> Initial state index for fast evolution
   integer, save :: FAST_INITIAL
   !> @}

   !> @name "Slow" Butcher Tableau
   !! @{
   integer, save :: ml_kmax
   real, dimension(:, :, :), allocatable, save :: ml_gamK, ml_wK
   real, dimension(:, :), allocatable, save :: ml_gamBar, ml_wBar
   real, dimension(:), allocatable, save :: ml_cS
   !> @}

   !> @name "Fast" Butcher Tableau
   !! @{
   real, dimension(:, :), allocatable, save :: ml_AF
   real, dimension(:), allocatable, save :: ml_bF, ml_cF
   !> @}

end module ml_mrData
