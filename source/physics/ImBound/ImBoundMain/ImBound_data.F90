!!****if* source/physics/ImBound/ImBoundMain/ImBound_data
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
!!  ImBound_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE ImBound_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the ImBound module.
!!
!!***
 
 
module ImBound_data

#include "Simulation.h"
#include "constants.h"

  logical, save :: ib_useImBound

  real, save :: ib_rhoSolid
  real, save :: ib_muSolid
  real, save :: ib_thcoSolid
  real, save :: ib_CpSolid

  integer, save :: ib_lsIt

  integer, save :: ib_meshMe
  integer, save :: ib_meshNumProcs
  integer, save :: ib_meshComm
  integer, save :: ib_nstep

  real, dimension(2), save :: ib_alfa
  real, save :: ib_rhoa
  real, save :: ib_gama

  integer, save :: ib_iVelFVar

end module ImBound_data
