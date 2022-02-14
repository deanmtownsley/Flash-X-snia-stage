!!****if* source/Particles/localAPI/pt_updateTypeDS
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
!!  pt_updateTypeDS
!!
!! SYNOPSIS
!!
!!  pt_updateTypeDS( integer(in)  :: particlesPerBlk(MAXBLOCKS,NPART_TYPES))
!!
!! DESCRIPTION
!!
!!  Part of the time advancement routine for the Particle Unit.
!!    Calculates the number of each type of particles in the sorted data structure 
!!    (optional return).    Returns the summed total of passive and active particles.
!!    Within the particles data structure, 
!!    PASSIVE_PART_TYPE in located in positions 1:totalPassive
!!    All active particle types are located in positions 
!!             totalPassive+1:totalPassive+1+totalActive
!!
!! ARGUMENTS
!!
!!  particlesPerBlk -- an array of size MAXBLOCKS,NPART_TYPES.  
!!                      It is produced by
!!                      Grid_sortParticles and holds the number 
!!                      of each type of particle in each block.
!!
!! NOTE
!!
!!   The values are actually calculated already in Grid_sortParticles but
!!   encapsulation rules make recalculating there here easier than passing them
!!   back and forth in data, etc.
!!
!!
!!***

!===============================================================================

subroutine pt_updateTypeDS (particlesPerBlk)

  implicit none

#include "Simulation.h"
  
  integer, INTENT(in), dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk


  return
  
end subroutine pt_updateTypeDS

!===============================================================================

