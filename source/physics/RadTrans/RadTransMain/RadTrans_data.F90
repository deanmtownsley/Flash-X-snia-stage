!!****if* source/physics/RadTrans/RadTransMain/RadTrans_data
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
!!
!!  RadTrans_data
!!
!!  SYNOPSIS
!!   use RadTrans_data
!!
!!  DESCRIPTION 
!!    Stores local data for the RadTrans unit
!!
!!***
module RadTrans_data

#include "constants.h"
#include "Simulation.h"

  implicit none
  
  integer, save :: rt_meshMe ! Process rank
  integer, save :: rt_meshNumProcs
  logical, save :: rt_useRadTrans
  logical, save :: rt_enableTiling

  integer, save :: rt_eosModeGc
  logical,dimension(NUNK_VARS),save :: rt_gcMask

  character(len=MAX_STRING_LENGTH) :: rt_str_geometry
  integer, save :: rt_geometry

  real, save :: rt_cfl
  
end module RadTrans_data
