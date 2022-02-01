!!****if* source/Particles/ParticlesMain/Particles_getCountPerBlk
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
!!  Particles_getCountPerBlk
!!
!! SYNOPSIS
!!
!!  Particles_getCountPerBlk(integer(OUT)  :: perBlkCount(:))
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Finds the number of particles in each block.  
!! 
!!
!! ARGUMENTS 
!!
!! perBlkCount : integer array containing number of particles on each blk.
!!
!!***

#include "Simulation.h"
subroutine Particles_getCountPerBlk(perBlkCount)

  use Particles_data, ONLY : particles, pt_numLocal
#ifdef DEBUG_PARTICLES
  use Logfile_interface, ONLY: Logfile_open, Logfile_close
#endif

  implicit none
#include "constants.h"

  integer,dimension(MAXBLOCKS),intent(OUT) :: perBlkCount
  integer :: i,j
#ifdef DEBUG_PARTICLES
  integer :: logunit
#endif

  perBlkCount=0
#ifdef DEBUG_PARTICLES
  call Logfile_open(logUnit, .TRUE.)
  write(logUnit,*) 'the number of particles',pt_numLocal
  call Logfile_close(.TRUE.)
#endif
  do i = 1,pt_numLocal
     j=int(particles(BLK_PART_PROP,i))
     if(j/=NONEXISTENT) then
        perBlkCount(j)=perBlkCount(j)+1
     end if
  end do
  
  return
end subroutine Particles_getCountPerBlk
