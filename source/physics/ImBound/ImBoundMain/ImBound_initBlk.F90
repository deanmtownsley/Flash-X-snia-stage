!!***if* source/physics/ImBound/ImBoundMain/ImBound_initBlk
!! NOTICE
!!  Copyright 2024 UChicago Argonne, LLC and contributors
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
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine ImBound_initBlk(tileDesc)

   use ImBound_data, ONLY: ib_numBodies
   use ImBound_type, ONLY: ImBound_type_t
   use ImBound_interface, ONLY: ImBound_getBodyPtr, ImBound_mapToGrid, ImBound_releaseBodyPtr
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Grid_tile, ONLY: Grid_tile_t
 
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   type(ImBound_type_t), pointer :: bodyInfo
   integer :: ibd
 
   nullify (bodyInfo)

   call Timers_start("ImBound_initBlk")

   do ibd = 1, ib_numBodies
      call ImBound_getBodyPtr(bodyInfo, ibd)
      call ImBound_mapToGrid(tileDesc, bodyInfo)
      call ImBound_releaseBodyPtr(bodyInfo, ibd)
   end do

   call Timers_stop("ImBound_initBlk")

end subroutine ImBound_initBlk
