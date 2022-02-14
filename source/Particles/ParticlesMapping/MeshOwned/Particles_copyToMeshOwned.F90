!!****if* source/Particles/ParticlesMapping/MeshOwned/Particles_copyToMeshOwned
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
!!  call Particles_copyToMeshOwned(particles2D)
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

subroutine Particles_copyToMeshOwned(particles2D, nprops, np, filePropNames)
    
  ! this is an alias for pt_amrexParticle_t in pt_amrexParticleMod:
  use Particles_dataTypes, ONLY : Particles_meshParticle_t
  use Particles_data, ONLY : pt_MeshMe, pt_containers, pt_flashProp2AoSRProp
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,        ONLY : Grid_tile_t

  use amrex_fort_module, only: amrex_long
  
  use pt_amrexParticleMod, ONLY : amrex_set_particle_id, amrex_set_particle_cpu

  implicit none

#include "constants.h"  
#include "Simulation.h"
#include "Particles.h"
  
    integer, INTENT(in) :: nprops, np
    real, dimension(1:nprops,1:np), INTENT(in) :: particles2D
    character(len=24), dimension(nprops), INTENT(in) :: filePropNames 
    type(Particles_meshParticle_t) :: meshOwnedPart
    type(Grid_tile_t)    :: tileDesc
    type(Grid_iterator_t) :: itor

    integer :: i,particleType, p_count, count, ii, c1
    integer :: j, d, cpu, ind
    integer(amrex_long) id
    integer :: grd_index,tile_index,level_amrex
!!------------------------------------------------------------------------------
   !Block information used by ParticleContainer functions
   call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
   if (itor%isValid()) then
      call itor%currentTile(tileDesc)
   
      level_amrex = tileDesc%level-1
      grd_index=tileDesc%grid_index
      tile_index=tileDesc%tile_index
 
      do j = 1, np
         ind = PART_TYPES_BEGIN
         do i = 1, nprops
            select case (filePropNames(i))
               case('posx')
                  meshOwnedPart%pos(1) = particles2D(i,j)
#if NDIM >= 2
               case('posy')
                  meshOwnedPart%pos(2) = particles2D(i,j)
#endif
#if NDIM == 3
               case('posz')
                  meshOwnedPart%pos(3) = particles2D(i,j)
#endif
               case('velx')
                  meshOwnedPart%vel(1) = particles2D(i,j)
#if NDIM >= 2
               case('vely')
                  meshOwnedPart%vel(2) = particles2D(i,j)
#endif
#if NDIM == 3
               case('velz')
                  meshOwnedPart%vel(3) = particles2D(i,j)
#endif
               case('cpu')
                  call amrex_set_particle_cpu(INT(particles2D(i,j)), meshOwnedPart)
               case('tag')
                  call amrex_set_particle_id(int(particles2D(i,j),kind=amrex_long), meshOwnedPart)
               case('proc')
                  ! do nothing
               case('blk')
                  ! do nothing 
               case('type')
                  particleType = int(particles2D(i,j))
                  if ((particleType >= PART_TYPES_BEGIN) .AND. &
                        (particleType <= PART_TYPES_END)) then
                     ind = particleType ! Change the slot in pt_containers where this particle goes
                  else
                     call Driver_abort("Particles_copyToMeshOwned: encountered invalid particle type")
                  end if
               case default
                  ! we should get here for pdens, ptemp etc.
                  if (pt_flashProp2AoSRProp(i) > 0 ) then
                     meshOwnedPart%fdata(pt_flashProp2AoSRProp(i)) = particles2D(i, j)
                  endif
            end select
         end do

         ! after setting all the properties for meshOwnedPart add it
         call pt_containers(ind)%add_particle(level_amrex, grd_index, tile_index, meshOwnedPart)
      end do
   endif
   call Grid_releaseTileIterator(itor)
   return

end subroutine Particles_copyToMeshOwned


