!!****if* source/Grid/GridParticles/GridParticlesMove/MeshOwned/Grid_countParticlesByBlock
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
!!  Grid_countParticlesByBlock
!!
!! SYNOPSIS
!!
!!  call Grid_countParticlesByBlock(particlesPerBlk)
!!                    
!!  
!! DESCRIPTION 
!!  
!!  counts the particles per block.
!! 
!!
!! ARGUMENTS 
!!
!!  particlesPerBlk:  Number of particles per block                      
!!
!! NOTES
!!   Currently this routine is called by io_ptwriteParticles (serial/parallel)
!!***

#ifdef DEBUG_ALL
#define DEBUG_SORT_PARTICLES
#endif

subroutine Grid_countParticlesByBlock(particlesPerBlk)

  use Driver_interface, ONLY : Driver_abort
  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_getTileIterator, Grid_releaseTileIterator, Grid_mapMeshToParticles

  use Particles_data, ONLY : pt_containers
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,        ONLY : Grid_tile_t

implicit none
#include "Simulation.h"
#include "constants.h"
#include "Particles.h"

  integer, dimension(MAXBLOCKS, NPART_TYPES) , intent(OUT)   :: particlesPerBlk
  integer :: p_count, ind, t_blk
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)    :: tileDesc
       
  ! Iterating thru blocks in exactly the same order that is used for writing UNK blocks to file
   call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
   t_blk = 1
   do while(itor%isValid())
        call itor%currentTile(tileDesc)
        do ind=1, NPART_TYPES
             p_count = pt_containers(ind)%num_particles(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index)
             ParticlesPerBlk(t_blk, ind)=p_count
             t_blk = t_blk+1
        enddo
        call itor%next()
   enddo             ! leaf itor enddo
   call Grid_releaseTileIterator(itor)
   return

end subroutine Grid_countParticlesByBlock
