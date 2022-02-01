!!****if* source/Grid/GridParticles/GridParticlesMove/MeshOwned/Grid_countParticles
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
!!  Grid_countParticles
!!
!! SYNOPSIS
!!
!!  call Grid_countParticles(integer(IN)    :: props
!!                          integer(INOUT) :: localCount,
!!                          integer(IN)    :: elementTypes,
!!                 optional,integer(IN)    :: maxCount)
!!                    
!!  
!! DESCRIPTION 
!!  
!!  counts the particles per type. 
!! 
!!
!! ARGUMENTS 
!!
!!  props : number of properties of each element in the dataBuf datastructure
!!
!!  localCount : While coming in it contains the current number of elements mapped to
!!                      this processor. After all the data structure movement, the number
!!                      of local elements might change, and the new value is put back into it.
!!  elementTypes  : Count of different types of elements in the simulation
!!  maxCount : This is parameter determined at runtime, and is the maximum number of local
!!                        elements that a simulation expects to have. All the arrays that hold
!!                        particles in the Particles unit are allocated based on this number.
!!                        
!!
!! NOTES
!!   Currently this routine is called by io_ptwriteParticles (serial/parallel), Grid_moveParticles and
!!   Particles_initData 
!!***

#ifdef DEBUG_ALL
#define DEBUG_SORT_PARTICLES
#endif

subroutine Grid_countParticles(props, localCount,elementTypes,maxCount)

  use Driver_interface, ONLY : Driver_abort
  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_getTileIterator, Grid_releaseTileIterator, Grid_mapMeshToParticles

  use Particles_data, ONLY : pt_containers
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,        ONLY : Grid_tile_t

implicit none
#include "Simulation.h"
#include "constants.h"
#include "Particles.h"

  integer, optional, intent(IN) ::  maxCount
  integer,intent(IN) ::  props,elementTypes
  integer, intent(INOUT) :: localCount
  integer,dimension(4, elementTypes)::localElementsPerType

  integer :: i, j,k,m,n,attrib,localNumBlocks,validParticles
  logical :: pType
  integer :: p_count, ind

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)    :: tileDesc

  localElementsPerType(:,:)=0

   call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  
  do while(itor%isValid())
      call itor%currentTile(tileDesc)
      do ind=1,elementTypes
        p_count = pt_containers(ind)%num_particles(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index)
        localElementsPerType(1,ind)=localElementsPerType(1,ind) + p_count
      enddo
      call itor%next()
   enddo             ! leaf itor enddo
   
   call Grid_releaseTileIterator(itor)

   localCount=sum(localElementsPerType(1,1:elementTypes))
   return

end subroutine Grid_countParticles
