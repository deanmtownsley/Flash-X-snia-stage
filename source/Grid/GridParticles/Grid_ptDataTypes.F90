#include "Simulation.h"
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

module Grid_ptDataTypes

#ifdef FLASH_GRID_AMREX

  use Particles_dataTypes, ONLY : Grid_particle_t => Particles_meshParticle_t
  use pt_amrexParticleMod, ONLY : Grid_ptContainer_t => pt_amrexParticlecontainer_t

#else

#  ifdef FLASH_USE_FORTRAN2003
  ! The Grid_ptContainer_t is currently not used anywhere. This is more
  ! for demonstration and experimentation with some newer Fortran features
  ! (specifically, parameterized derived types) than anything else.
  type Grid_ptContainer_t (nParticles)
     integer, LEN :: nParticles = 1
     real,dimension(1:NPART_PROPS,nParticles) :: data
  end type Grid_ptContainer_t
#  else
  type Grid_ptContainer_t
     real,dimension(:,:),allocatable :: data
  end type Grid_ptContainer_t
#  endif

  type Grid_particle_t
     real,dimension(1:NPART_PROPS) :: data
  end type Grid_particle_t

#endif
end module Grid_ptDataTypes
