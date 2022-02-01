!!****if* source/Particles/ParticlesMapping/MeshOwned/Particles_copyFromMeshOwned
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
!!  Grid_copyFromMeshOwned
!!
!! SYNOPSIS
!!
!!  call Particles_copyFromMeshOwned(np, particleType, meshOwnedPart, particles2D)
!!                    
!!  
!! DESCRIPTION 
!!  
!! Get position, velocity, id, cpu from amrex particle and 
!! populate particles2D for writing into hdf5 file later
!!
!! ARGUMENTS 
!!
!!  np : number of particles
!!
!!  meshOwnedPart : particles in mesh-owned format
!!
!!  particles2D : 2D particles with property values
!!
!!  particleType : Type of particle : passive, active, star etc.
!!
!!***
!===============================================================================

subroutine Particles_copyFromMeshOwned(np, particleType, meshOwnedPart, particles2D)
    
  ! this is an alias for amrex_particle in pt_amrexParticleMod:
  use Particles_dataTypes, ONLY : Particles_meshParticle_t

  use Particles_data, ONLY : pt_MeshMe, pt_numMappedAttributes, pt_mappedAttributes, pt_aoSRProp2FlashProp

  use pt_amrexParticleMod, ONLY : amrex_get_particle_cpu, amrex_get_particle_id
  use amrex_fort_module, only: amrex_long
  
  implicit none

#include "constants.h"  
#include "Simulation.h"
#include "Particles.h"
  
  integer, INTENT(in) :: np, particleType
  type(Particles_meshParticle_t), INTENT(in) :: meshOwnedPart(:)
  real, dimension(:,:), INTENT(out)          :: particles2D

  integer             :: i,particleTypes, p_count, count
  integer             :: j, d, cpu, k
  integer(amrex_long) :: id
!!------------------------------------------------------------------------------
   do j=1, np
      ! set default value to -1
      particles2D(:,j) = -1
      do d=1, NDIM
         particles2D(POSX_PART_PROP+d-1, j) = meshOwnedPart(j)%pos(d)
         particles2D(VELX_PART_PROP+d-1, j) = meshOwnedPart(j)%vel(d)
      enddo
      ! set remaining pos/vel to zero, if ndim is 1 set 2/3 to 0.
      ! if ndim is 2 set 3 to 0.
      if (NDIM .lt. 3) then
         do d=NDIM+1, 3
            particles2D(POSX_PART_PROP+d-1, j) = 0.0
            particles2D(VELX_PART_PROP+d-1, j) = 0.0
         enddo
      endif

#ifdef TYPE_PART_PROP
  particles2D(TYPE_PART_PROP, j) = particleType
#endif

      ! set tag as the id, cpu from amrex as cpu (birthplace of particle) and 
      ! proc as the current processor of the particle
      call amrex_get_particle_cpu(cpu, meshOwnedPart(j))
      call amrex_get_particle_id(id, meshOwnedPart(j))

      particles2D(TAG_PART_PROP, j) = REAL(id)
      particles2D(CPU_PART_PROP, j) = REAL(cpu)
      particles2D(PROC_PART_PROP, j) = REAL(pt_MeshMe)
   enddo

   ! moved to a separate loop going from k->j for better performance
   ! handle simulation specific attributes such as ptemp/pdens 
   do k=1, size(pt_aoSRProp2FlashProp)  
      do j=1, np
         ! print *, pt_MeshMe, j, k, aoSRProp2FlashProp(k)
         if (pt_aoSRProp2FlashProp(k) > 0 ) then ! check for index > 0 
            particles2D(pt_aoSRProp2FlashProp(k), j) = meshOwnedPart(j)%fdata(k)
         endif
      enddo
   enddo

return

end subroutine Particles_copyFromMeshOwned


