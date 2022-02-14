!!****if* source/Particles/ParticlesMain/pt_initFlashAoSRMap
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
!!  pt_initFlashAoSRMap
!!
!! SYNOPSIS
!!
!!  call pt_initFlashAoSRMap()
!!
!! DESCRIPTION
!!
!!  Sets the pt_flashProp2AoSRProp and pt_aoSRProp2FlashProp mapping.
!!
!!
!! NOTES
!! used in Particles_copyToMeshOwned
!!***
!===============================================================================

subroutine pt_initFlashAoSRMap ()

   use Particles_data, ONLY: pt_flashProp2AoSRProp, pt_aoSRProp2FlashProp
   use Simulation_interface, ONLY : Simulation_mapIntToStr
  
   implicit none
 
#include "constants.h"
#include "Simulation.h"
#include "Particles.h"
   
   integer :: i, count
   count = 1
   ! initialize default to -1
   pt_aoSRProp2FlashProp(:) = -1
   pt_flashProp2AoSRProp(:) = -1

   do i=1, NPART_PROPS
#if defined(CPU_PART_PROP)
      if(i == CPU_PART_PROP) then 
         cycle
      endif
#endif
#if defined(POSX_PART_PROP) || defined(POSY_PART_PROP) || defined(POSZ_PART_PROP)
      if((i == POSX_PART_PROP) .or. (i == POSY_PART_PROP) .or. (i == POSZ_PART_PROP)) then 
         cycle
      endif
#endif
#if defined(VELX_PART_PROP) || defined(VELY_PART_PROP) || defined(VELZ_PART_PROP)
      if((i == VELX_PART_PROP) .or. (i == VELY_PART_PROP) .or. (i == VELZ_PART_PROP)) then 
         cycle
      endif
#endif
#if defined(TYPE_PART_PROP)
      if(i == TYPE_PART_PROP) then 
         cycle
      endif
#endif
#if defined(BLK_PART_PROP)
      if(i == BLK_PART_PROP) then 
         cycle
      endif
#endif
#if defined(PROC_PART_PROP)
      if(i == PROC_PART_PROP) then 
         cycle
      endif
#endif
#if defined(TAG_PART_PROP)
      if(i == TAG_PART_PROP) then 
         cycle
      endif
#endif

#ifdef DEBUG_PARTICLES
         print *, '[pt_initFlashAoSRMap]:', i,' FlashPProp to AoSRProp ', count
         print *, '[pt_initFlashAoSRMap]:', count, ' AoSRProp to FlashPProp ', i 
#endif         
      pt_flashProp2AoSRProp(i) = count
      pt_aoSRProp2FlashProp(count) = i
      count = count + 1
   enddo
   
 end subroutine pt_initFlashAoSRMap
 
 
 
