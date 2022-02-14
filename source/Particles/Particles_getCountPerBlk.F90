!!****f* source/Particles/Particles_getCountPerBlk
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
subroutine Particles_getCountPerBlk(perBlkCount)

  implicit none
#include "Simulation.h"

  integer,dimension(MAXBLOCKS),intent(OUT) :: perBlkCount

  perBlkCount=0
  
  return
end subroutine Particles_getCountPerBlk
