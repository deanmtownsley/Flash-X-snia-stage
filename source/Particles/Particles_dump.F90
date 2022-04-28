!!****f* source/Particles/Particles_dump
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
!!    Particles_dump
!!
!! SYNOPSIS
!!    Particles_dump()
!!                   integer(IN) :: blockCount,
!!                   integer(IN), dimension(:)  :: blockList,
!!                   integer(IN) :: nstep,
!!                   real(IN)    :: time)
!!
!! DESCRIPTION
!!
!!   Dump particle information to a plain old file.  The output file is called
!!      test_ParticlesDump_00##  where the ## relates to the processor number.
!!   Quick and dirty output implemented for StirTurb BG/L testing.  You can see an
!!      example of the output by running the simulation unitTest/Particles.
!!   On the positive side, there is an associated fidlr3/IDL routine particles_dump.pro which will
!!      read the input from this file.
!!
!! ARGUMENTS
!!
!!  
!!  blockCount:             integer(IN)     number of blocks within this processor
!!  blockList(blockCount):  integer(IN)     block IDs within this processor
!!  nstep:                  integer(IN)     current time step index
!!  time:                   real(IN)        current simulation time
!!
!! PARAMETERS
!!
!! NOTES
!!
!!  The source for the real routine is in the ParticlesMain/unitTest subdirectory
!!
!!***


!-------------------------------------------------------------------
subroutine Particles_dump(blockCount,blockList,nstep,time,dt)

#include "Simulation.h"
#include "constants.h"

  implicit none

  
  integer, intent(IN) :: blockCount
  integer, intent(IN) :: blockList(blockCount)
  integer, intent(IN) :: nstep
  real, intent(IN)    :: time, dt
 
  
end subroutine Particles_dump

