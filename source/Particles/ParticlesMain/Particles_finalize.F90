!!****if* source/Particles/ParticlesMain/Particles_finalize
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
!!    Particles_finalize
!!
!! SYNOPSIS
!!    Particles_finalize( )
!!
!! DESCRIPTION
!!
!!    Finalize routine for the particle unit.  Removes memory usage
!!      set in Particles_init.
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!
!!***

#include "Simulation.h"

subroutine Particles_finalize()
  use pt_interface, ONLY : pt_initFinalize
  use Particles_data, ONLY : particles,useParticles

#ifdef FLASH_GRID_AMREX
  use Particles_data, ONLY : pt_containers
  use pt_amrexParticleMod, ONLY: amrex_particlecontainer_destroy
#endif
  implicit none


  integer :: i
  if(useParticles) deallocate(particles)
#ifdef FLASH_GRID_AMREX
  if(useParticles) then
    do i=1,NPART_TYPES
       call amrex_particlecontainer_destroy(pt_containers(i))
    end do
  end if
#endif
  call pt_initFinalize()
end subroutine Particles_finalize
