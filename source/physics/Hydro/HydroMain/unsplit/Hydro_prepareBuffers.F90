!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_prepareBuffers
!!
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
!!
!! NAME
!!
!!  Hydro_prepareBuffers
!!  For more details see the documentation of the NULL implementation
!!
!!***



#include "Simulation.h"
subroutine Hydro_prepareBuffers
  use hy_memInterface, ONLY :  hy_memAllocScratch
  use Hydro_data, ONLY : hy_fullRiemannStateArrays
  implicit none

#include "constants.h"
#include "UHD.h"

#ifndef FLASH_GRID_AMREX
#ifndef FLASH_GRID_UG
  integer :: blockList(MAXBLOCKS)
  integer :: blockCount
#endif


  ! Independently of whether hy_fullRiemannStateArrays is set:
  call hy_memAllocScratch(SCRATCH_CTR,HY_VAR1_SCRATCHCTR_VAR,2, 0,0,0)
!!$       blockList(1:blockCount) )

#ifndef FLASH_GRID_UG
  if (hy_fullRiemannStateArrays) then
     call hy_memAllocScratch(SCRATCH_FACEX,&
       min(HY_P01_FACEXPTR_VAR,HY_N01_FACEXPTR_VAR),&
       2*HY_SCRATCH_NUM, &
       0,1,0)
!!$     , &
!!$       blockList(1:blockCount) )
#if NDIM > 1
     call hy_memAllocScratch(SCRATCH_FACEY,&
       min(HY_P01_FACEYPTR_VAR,HY_N01_FACEYPTR_VAR),&
       2*HY_SCRATCH_NUM,&
       0,1,0)
!!$     , &
!!$       blockList(1:blockCount) )
#endif
#if NDIM > 2
     call hy_memAllocScratch(SCRATCH_FACEZ,&
       min(HY_P01_FACEZPTR_VAR,HY_N01_FACEZPTR_VAR),&
       2*HY_SCRATCH_NUM,&
       0,1,0)
!!$     , &
!!$       blockList(1:blockCount) )
#endif
  end if

#endif
#endif
  
end subroutine Hydro_prepareBuffers
