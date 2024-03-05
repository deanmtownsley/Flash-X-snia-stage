!!****if* source/Grid/GridParticles/GridParticlesMapFromMesh/Grid_mapMeshToParticles
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
!!  Grid_mapMeshToParticles(real(inout) :: particles(part_props,numParticles),
!!                          integer(in) :: part_props,
!!                          integer(in) :: part_blkID,
!!                          integer(in) :: numParticles,
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
!!     particles:     Data structure containing particles information
!!     part_props : number of particle attributes
!!     part_blkID : the index of particle attributes that carries the block number
!!     numParticles : the number of particles on my proc
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
!!
!!  
!!***
!*******************************************************************************
#include "Simulation.h"
#include "constants.h"
#include "GridParticles.h"
#include "Particles.h"

subroutine Grid_mapMeshToParticles (particles, part_props,part_blkID,&
                                    numParticles,posAttrib,&
                                    numAttrib, attrib,&
                                    mapType,gridDataStruct)

  use Driver_interface, ONLY : Driver_abort
  use Particles_interface, ONLY : Particles_mapFromMesh
  use physicaldata, ONLY : unk
  use Grid_data, ONLY : gr_delta
#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : bnd_box, lrefine, lnblocks
#endif
  implicit none



  integer, INTENT(in) :: part_props, numParticles, numAttrib, part_blkID
  real, INTENT(inout),dimension(part_props,numParticles) :: particles
  integer,dimension(MDIM), intent(IN) :: posAttrib
  integer, dimension(PART_ATTR_DS_SIZE,numAttrib),INTENT(in) :: attrib
  integer, INTENT(IN) :: mapType
  integer, optional, intent(IN) :: gridDataStruct

  integer :: i,j,k,blk

  real, pointer, dimension(:,:,:,:) :: solnVec
  real, dimension(LOW:HIGH,MDIM) :: bndBox
  real,dimension(MDIM) :: delta, pos
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits
  integer :: gDataStruct
  real, dimension(numAttrib) :: partAttribVec
  integer :: level, blkCount
  

  
  if(present(gridDataStruct)) then
     gDataStruct=gridDataStruct
  else
     gDataStruct=CENTER
  end if

  
  nullify(solnVec)

  if(numParticles>0) then

     k=0
     blkLimits=1
     blkLimits(LOW,1:NDIM)=1+NGUARD
     blkLimits(HIGH,IAXIS)=blkLimits(LOW,IAXIS)+NXB
     blkLimits(HIGH,IAXIS)=blkLimits(LOW,JAXIS)+NYB*K2D
     blkLimits(HIGH,IAXIS)=blkLimits(LOW,KAXIS)+NZB*K3D
     
     do i=1,numParticles 

#ifdef DEBUG_GRIDPARTICLES
        if((particles(part_blkID, i) < 0) .or. (particles(part_blkID, i) > blkCount)) then
           call Driver_abort("BLK_PART_PROP out of bounds")
        end if
#endif
#ifdef FLASH_GRID_PARAMESH
        blk = int(particles(part_blkID,i))
        blkCount=lnblocks
        level=lrefine(blk)
        bndBox(LOW:HIGH,1:MDIM)=bnd_box(LOW:HIGH,1:MDIM,blk)
#else
        blkCount=1
        blk=1
        level=1
        call gr_getBndBox(bndBox)
#endif        
        solnVec(1:,1:,1:,1:)=>unk(:,:,:,:,blk)
        delta(:)=gr_delta(:,level)
        pos(1:MDIM)=particles(posAttrib(1:MDIM),i)
        call Particles_mapFromMesh (mapType, numAttrib, attrib,&
             pos, bndBox,delta,blkLimits,solnVec, partAttribVec)
        particles(attrib(PART_DS_IND,1:numAttrib),i)=partAttribVec(1:numAttrib)
     end do
  end if
  
end subroutine Grid_mapMeshToParticles

