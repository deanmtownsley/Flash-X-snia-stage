!!****if* source/numericalTools/MoL/MoLMain/ERK/ml_data
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
!!    Stores data for the explicit Runge-Kutta (ERK) integrator
!!
!!  NOTES
!!
!!***
module ml_data

   implicit none

   character(len=:), allocatable, save :: ml_method

   integer, save :: ml_order
   integer, save :: ml_stages

   ! RK tableau
   real, dimension(:, :), allocatable, target, save :: ml_A
   real, dimension(:), allocatable, target, save :: ml_b
   real, dimension(:), allocatable, target, save :: ml_c

   ! Indices for intermediate RHS states
   integer, allocatable, save :: ml_K(:)

end module ml_data
