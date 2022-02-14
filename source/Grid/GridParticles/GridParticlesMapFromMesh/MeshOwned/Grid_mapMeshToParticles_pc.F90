!!****f* source/Grid/GridParticles/GridParticlesMapFromMesh/Amrex
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
!!  Grid_mapMeshToParticles
!!
!! SYNOPSIS
!!
!!  Grid_mapMeshToParticles(int(in) :: ptContainerPos,
!!                          integer(in) :: part_props,
!!                          integer(in) :: part_blkID,
!!                          integer(in) :: posAttrib(POS),
!!                          integer(in) :: numAttrib,
!!                          integer(in) :: attrib(:,numAttrib),
!!                          integer(in) :: mapType,
!!               optional,  integer(in) :: gridDataStruct)
!!
!! DESCRIPTION
!!
!!  Routine to map a quantity defined on the mesh onto the particle positions.
!!
!! ARGUMENTS
!!
!!     ptContainerPos:     Position of particles type in AMReX's ParticleContainer
!!                                    array pt_containers of Particle_data module
!!     part_props : number of particle attributes
!!     part_blkID : the index of particle attributes that carries the block number (only for paramesh)
!!     posAttrib           : particles data structure indices pointing
!!                           to particle positions
!!     numAttrib           : number of attributes that need to be mapped
!!     attrib              : list of attributes and their corresponding
!!                           mesh data structure indices
!!                           processor
!!     mapType : method for mapping grid quantities to particles
!!     gridDataStruct: The Grid data structure that varGrid refers to; one of
!!                    CENTER, FACEX, FACEY, FACEZ, GRIDVAR.  
!!                    If this argument is not present, the default is CENTER, 
!!                    so that attrib indicates a list of variable in UNK.
!! NOTES: 
!! Arguments part_blkID is dummy for amrex implementation
!!  
!!***
!*******************************************************************************

#include "Simulation.h"
#include "constants.h"
#include "Particles.h"

subroutine Grid_mapMeshToParticles_pc (ptContainerPos, part_props,part_blkID,&
                                    posAttrib,&
                                    numAttrib, attrib,&
                                    mapType,gridDataStruct)

  use Driver_interface, ONLY : Driver_abort
  ! this is an alias for pt_amrexParticle_t in pt_amrexParticleMod:
  use Grid_ptDataTypes, ONLY : Grid_particle_t
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator
  use Particles_interface, ONLY : Particles_mapFromMesh
  use Particles_data, ONLY : pt_containers
#ifdef DEBUG_PARTICLES
  use Particles_data, ONLY : pt_meshMe
#endif

  use pt_interface, ONLY : pt_initFlashAoSRMap
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,        ONLY : Grid_tile_t
  use Particles_data, ONLY: pt_flashProp2AoSRProp

  implicit none

  integer, INTENT(in) :: part_props, numAttrib, part_blkID
  integer, INTENT(IN) :: ptContainerPos
  integer,dimension(MDIM), intent(IN) :: posAttrib
  integer, dimension(PART_ATTR_DS_SIZE,numAttrib),INTENT(in) :: attrib
  integer, INTENT(IN) :: mapType
  integer, optional, intent(IN) :: gridDataStruct
  !-----------------------------------------------------------------------------------------------!
  integer :: gDataStruct, count
  real, contiguous, pointer :: solnData(:,:,:,:)
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)    :: tileDesc
  type(Grid_particle_t), pointer :: particles(:)
  integer :: numParticlesOnBlock, i, j, prop_id, jAoSR
#ifdef DEBUG_PARTICLES
  logical :: map_velocity
#endif
  real, dimension(LOW:HIGH,MDIM) :: bndBox
  real,dimension(MDIM) :: delta, pos
  real, dimension(numAttrib) :: partAttribVec

  nullify(solnData)
  if(present(gridDataStruct)) then
     gDataStruct=gridDataStruct
  else
     gDataStruct=CENTER
  end if
!! Using grid iterator will return all blocks/tiles.  Even the ones that do not have 
!! particles on this container. Better way is to use the "ParIter'' which is not
!! yet exposed from AMReX in Fortran. There is no wrapper in FLASH for the same.
  count = 0

#ifdef DEBUG_PARTICLES
  ! the index attrib(PART_DS_IND,numAttrib) refers to the last index for mapping
  ! checking to see if that index is for AoSRProp or fdata properties
  ! if not, we'll assume this is processing velocities. 
  prop_id = attrib(PART_DS_IND,numAttrib)
  map_velocity = .true.
 
  if (pt_flashProp2AoSRProp(prop_id) .ne. -1) then
     map_velocity = .false.
  endif

  print *, "Debug [Grid_mapMeshToParticles_pc]:", pt_meshMe, "prop_id, map_velocity?", prop_id, map_velocity
#endif                          

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)

   do while(itor%isValid())
      call itor%currentTile(tileDesc)
      call tileDesc%getDataPtr(solnData, CENTER)
!!TODO:: Particles_getFromTile(pt_containers(ptContainerPos), tileDesc, particles)   !!Correct usage once a wrapper get function is introduced
      particles => pt_containers(ptContainerPos)%get_particles(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index)
      numParticlesOnBlock = size(particles)
      count = count + numParticlesOnBlock

      if(numParticlesOnBlock>0) then
         call tileDesc%boundBox(bndBox)
         call tileDesc%deltas(delta)
         do i = 1, numParticlesOnBlock
            do j = 1,NDIM
               pos(j) = particles(i)%pos(j)
            end do
            call Particles_mapFromMesh (mapType, numAttrib, attrib,&
            pos, bndBox,delta, tileDesc%limits, solnData, partAttribVec)

            do j=1, numAttrib
               prop_id = attrib(PART_DS_IND,j)
               if (prop_id > 0) then
                  select case (prop_id)
                  case(POSX_PART_PROP)
                     particles(i)%pos(1) = partAttribVec(j)
#if NDIM >= 2
                  case(POSY_PART_PROP)
                     particles(i)%pos(2) = partAttribVec(j)
#endif
#if NDIM == 3
                  case(POSZ_PART_PROP)
                     particles(i)%pos(3) = partAttribVec(j)
#endif
                  case(VELX_PART_PROP)
                     particles(i)%vel(1) = partAttribVec(j)
#if NDIM >= 2
                  case(VELY_PART_PROP)
                     particles(i)%vel(2) = partAttribVec(j)
#endif
#if NDIM == 3
                  case(VELZ_PART_PROP)
                     particles(i)%vel(3) = partAttribVec(j)
#endif
                  case default
                     jAoSR = pt_flashProp2AoSRProp(prop_id)
                     if (jAoSR > 0) particles(i)%fdata(jAoSR) = partAttribVec(j)
                     ! Maybe add error handling if jAoSR<=0 ? - KW
                  end select
#ifdef DEBUG_PARTICLES
               else
                  call Driver_abort("Grid_mapMeshToParticles_pc: prop_id is invalid, aborting.")
#endif
               end if
            enddo

         end do
      end if

      call tileDesc%releaseDataPtr(solnData, CENTER)
      nullify(solnData)
      call itor%next()
   enddo
#ifdef DEBUG_PARTICLES
   print *, "Debug:", pt_meshMe, "numparticles mapped: ", count
#endif
  call Grid_releaseTileIterator(itor)

  return
end subroutine Grid_mapMeshToParticles_pc
