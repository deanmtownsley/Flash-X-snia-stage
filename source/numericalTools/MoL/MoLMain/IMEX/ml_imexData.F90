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
!! @brief Main data module for MoL's IMEX integrator

!> @ingroup MoLMR
!! Stores data for the IMEX integrator
module ml_imexData

   implicit none

   !> @name IMEX Method Information
   !! @{

   !> Name of the IMEX method
   character(len=:), allocatable, save :: ml_method

   !> Order of the IMEX method
   integer, save :: ml_order

   !> Number of stages
   integer, save :: ml_nstages

   !> @}

   !> @name Indexing
   !! @{
   !! @note Indexing variables below do not contain an ml_ to maintain a consistent notation
   !!       with the pre-processor defined variable indexing in Simulation.h and elsewhere

   !> Explicit RHS indices
   integer, dimension(:), allocatable, save :: FE

   !> Implicit RHS indices
   integer, dimension(:), allocatable, save :: FI

   !> @}

   !> @name Explicit Butcher Tableau
   !! @{
   real, dimension(:, :), allocatable, target, save :: ml_AE
   real, dimension(:), allocatable, save :: ml_bE, ml_cE
   !> @}

   !> @name Implicit Butcher Tableau
   !! @{
   real, dimension(:, :), allocatable, target, save :: ml_AI
   real, dimension(:), allocatable, save :: ml_bI, ml_cI
   !> @}

end module ml_imexData
