!!****if* source/numericalTools/MoL/MoLMain/MR/mr_data
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
!!    mr_data
!!
!!  SYNOPSIS
!!    use mr_data
!!
!!  DESCRIPTION
!!    Stores data for the multi-rate implicit-explicit integrator
!!
!!  NOTES
!!
!!***
module mr_data

   implicit none

   character(len=:), allocatable, save :: mr_slowMethod
   character(len=:), allocatable, save :: mr_fastMethod

   integer, save :: mr_nstages_slow, mr_nstages_fast
   integer, save :: mr_nsubcycle

   integer, save :: mr_nscratch_slow, mr_nscratch_fast

   ! RHS indexinng
   integer, dimension(:), allocatable, save :: FE, FI, FF

   ! Initial state for fast evolution
   integer, save :: FAST_INITIAL

   ! Butcher tableau for slow evolution
   integer, save :: kmax
   real, dimension(:, :, :), allocatable, save :: gamK, wK
   real, dimension(:, :), allocatable, save :: gamBar, wBar
   real, dimension(:), allocatable, save :: cS

   ! Butcher tableau for fast evolution
   real, dimension(:, :), allocatable, save :: AF
   real, dimension(:), allocatable, save :: bF, cF

end module mr_data
