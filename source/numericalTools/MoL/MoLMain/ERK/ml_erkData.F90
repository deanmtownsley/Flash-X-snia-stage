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
!! @brief Main data module for MoL's ERK integrator

!> @ingroup MoLERK
!! Stores data for the explicit Runge-Kutta (ERK) integrator
module ml_erkData

   implicit none

   !> @name ERK Method Information
   !! @{

   !> Name of the ERK method
   character(len=:), allocatable, save :: ml_method

   !> Order of the ERK method
   integer, save :: ml_order

   !> Number of stages in the ERK method
   integer, save :: ml_stages
   !> @}

   !> @name RK tableau
   !! @{

   !> Per-stage weights
   real, dimension(:, :), allocatable, target, save :: ml_A

   !> Weighting coefficients for final combination of intermedate statea
   real, dimension(:), allocatable, target, save :: ml_b

   !> Timing coefficients
   real, dimension(:), allocatable, target, save :: ml_c
   !> @}

   !> Indices for intermediate RHS states
   integer, allocatable, save :: ml_K(:)

end module ml_erkData
