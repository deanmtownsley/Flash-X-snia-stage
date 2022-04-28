!!****if* source/Grid/GridParticles/GridParticlesMove/paramesh/VirtualParticles/gr_ptVPData
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
!!  gr_ptVPData
!!
!! SYNOPSIS
!!
!!  use gr_ptVPData
!!
!! DESCRIPTION
!!
!!  Data module for the variables when using virtual particles
!!  with paramesh
!!
!! ARGUMENTS
!!
!!
!!***

#include "constants.h"
Module gr_ptVPData
  real,dimension(LOW:HIGH,MDIM),save :: gr_ptVPBndBox
  real,dimension(MDIM):: gr_ptVPDeltas
  integer, save :: gr_ptVPMaxCount=100
  real, save :: gr_ptVPBufferFactor=1.
end Module gr_ptVPData
