!!****if* source/IO/IOParticles/IOParticles_data
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
!!  IOParticles_data
!!
!! SYNOPSIS
!!
!!  IOParticles_data()
!!
!! DESCRIPTION 
!!  
!!  Holds all the IO data that is needed by the IOParticles unit  
!!
!!
!!***

module IOParticles_data

 
#include "constants.h"
#include "Simulation.h"


  integer, save :: io_particleFileNumber
  integer, save :: io_particleFileIntervalStep
  real, save :: io_particleFileIntervalTime
  integer, save :: io_nextParticleFileStep
  real, save :: io_nextParticleFileTime
  
  real, save :: io_particleFileIntervalZ
  real, save :: io_nextParticleFileZ

  logical, save :: useParticles
  logical, save :: io_dumpParticleFileExist
  logical, save :: io_ptSupressSinglePrecision
  logical, save :: io_writeParticleSubset
  logical, save :: io_writeParticleAll

end module IOParticles_data
