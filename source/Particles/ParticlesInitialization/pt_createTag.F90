!!****if* source/Particles/ParticlesInitialization/pt_createTag
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
!!    pt_createTag
!!
!! SYNOPSIS
!!
!!    call pt_createTag()
!!
!! DESCRIPTION
!!
!!    Initialize unique tag attribute values for all particle.
!!
!!    This algorithm works by first assigning tag values
!!    starting from one in each block. At the end of this step
!!    each LEAF block has a count of particles it is getting, and
!!    the tags are duplicated on each block. In the second step
!!    we add up all particles in the blocks to the left of my
!!    block and add that number to each particle's tag. At this
!!    point each particle in a processor has a unique tag, but
!!    the tag numbers are duplicated across processors. In the third
!!    and final step, using MPI_Allgather, the sum of the 
!!    number of particles in the processor left of pt_meshMe
!!    is added to current tag in each particle, thus generating
!!    unique tag for all particles globally.
!!
!! NOTES
!!    This method of tag generation will work for up to 10^14 
!!    particles in a simulation.
!!
!!
!!
!!
!!***

!!#define DEBUG_PARTICLES

subroutine pt_createTag()

  use Particles_data, ONLY : pt_numLocal,particles,pt_meshMe,pt_meshNumProcs,&
       pt_meshComm
  use Driver_interface, ONLY : Driver_abort
  use Grid_interface, ONLY : Grid_getLocalNumBlks
  use pt_interface, ONLY : pt_findTagOffset
  implicit none

#include "Simulation.h"
#include "constants.h"
  
  integer :: localNumBlocks,blockID,i, count, ierr, MyTagStart
  
  integer, allocatable ::  myBlkParticles(:) 
  
#ifndef PARTICLES_BIRTHCPU
  call Grid_getLocalNumBlks(localNumBlocks)

! #ifdef DEBUG_PARTICLES
!   if (pt_meshMe == MASTER_PE .OR. pt_meshNumProcs .LE. 4) then
!      print*,'pt_createTag on processor', pt_meshMe, ': localNumBlocks=',localNumBlocks
!   end if
! #endif

  allocate(myBlkParticles(localNumBlocks))

  myBlkParticles(1:localNumBlocks) = 0

  do i = 1,pt_numLocal
     blockID = int(particles(BLK_PART_PROP,i))
! #ifdef DEBUG_PARTICLES
!   if (pt_meshMe == MASTER_PE .OR. pt_meshNumProcs .LE. 4) then
!      print*,pt_meshMe, ': particle #', i,' on block #',blockID
!   end if
! #endif
     myBlkParticles(blockID)=myBlkParticles(blockID)+1
     particles(TAG_PART_PROP,i)=real(myBlkParticles(blockID))
  end do
  
  do i = 2,localNumBlocks
     myBlkParticles(i)=myBlkParticles(i)+myBlkParticles(i-1)
  end do
  

  call pt_findTagOffset(pt_numLocal,MyTagStart)

  do i=1,pt_numLocal
     blockID = int(particles(BLK_PART_PROP,i))
     if(blockID>1) then
        particles(TAG_PART_PROP,i)= particles(TAG_PART_PROP,i)+ &
             real(myBlkParticles(blockID-1)+&
             MyTagStart)
     else
        particles(TAG_PART_PROP,i)=particles(TAG_PART_PROP,i)+&
             real(MyTagStart)
     end if
  end do

  deallocate(myBlkParticles)
#else
   do i=1,pt_numLocal
      particles(TAG_PART_PROP,i)= i
      particles(CPU_PART_PROP,i)= pt_meshMe
   enddo
#endif
  return
end subroutine pt_createTag
