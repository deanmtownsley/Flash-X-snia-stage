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

subroutine Particles_createDataStructs
  use Particles_data, ONLY :  pt_containers, useParticles
  use amrex_amr_module, ONLY : amrex_get_amrcore
  use pt_amrexParticleMod, ONLY : amrex_particlecontainer_build

  implicit none
    integer :: i
    if(useParticles) then
        do i=1, NPART_TYPES
        call amrex_particlecontainer_build(pt_containers(i), amrex_get_amrcore())
        end do
    endif

end subroutine Particles_createDataStructs

