!!****if* source/numericalTools/MoL/MoLMain/ERK/erk_data
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
!!    erk_data
!!
!!  SYNOPSIS
!!    use erk_data
!!
!!  DESCRIPTION
!!    Stores data for the explicit Runge-Kutta (ERK) integrator
!!
!!  NOTES
!!
!!***
module erk_data

   implicit none

   character(len=:), allocatable, save :: erk_method

   integer, save :: erk_order
   integer, save :: erk_stages

   ! RK tableau
   real, dimension(:, :), allocatable, target, save :: erk_A
   real, dimension(:), allocatable, target, save :: erk_b
   real, dimension(:), allocatable, target, save :: erk_c

   ! Indices for intermediate RHS states
   integer, allocatable, save :: erk_K(:)

end module erk_data
