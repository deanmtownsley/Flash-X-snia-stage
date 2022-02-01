!!****f* source/Grid/Amrex
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
!!  Grid_mapMeshToParticles(int(in) :: pt_containerPos,
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
!!     pt_containerPos:     Position of particles type in AMReX's ParticleContainer
!!                                    array pt_containers of Particle_data module
!!     part_props : number of particle attributes
!!     part_blkID : the index of particle attributes that carries the block number
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

subroutine Grid_mapMeshToParticles_pc (pt_containerPos, part_props,part_blkID,&
                                    posAttrib,&
                                    numAttrib, attrib,&
                                    mapType,gridDataStruct)

#include "constants.h"
#include "Particles.h"
  implicit none

  integer, INTENT(in) :: part_props, numAttrib, part_blkID
  integer, INTENT(IN) :: pt_containerPos
  integer,dimension(MDIM), intent(IN) :: posAttrib
  integer, dimension(PART_ATTR_DS_SIZE,numAttrib),INTENT(in) :: attrib
  integer, INTENT(IN) :: mapType
  integer, optional, intent(IN) :: gridDataStruct

  return
end subroutine Grid_mapMeshToParticles_pc
