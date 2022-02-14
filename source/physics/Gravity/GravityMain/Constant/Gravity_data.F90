!!****if* source/physics/Gravity/GravityMain/Constant/Gravity_data
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
!!
!!  Gravity_data
!!
!! SYNOPSIS
!!
!!  use Gravity_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for constant gravity.
!!
!! PARAMTERS
!!
!!   grv_direc , string : specifies the direction of constant gravity
!!   grv_vector : value of constant gravity
!!
!!***

module Gravity_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!

  real, save :: grv_const
  character (len=80), save :: grv_direc

  !! *** Module Variables *** !!

  real, dimension(MDIM), save :: grv_vector
  integer, save :: grv_meshMe, grv_meshNumProcs
  logical, save :: useGravity

end module Gravity_data
