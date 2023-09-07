!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! 
!! 
!! This is AMReX specific implementation. For now it is identical to the
!! the one in GridMain, but in future if we allow use of non octree use of
!! AMReX it will have to be modified.
!!


#include "Simulation.h"
subroutine Grid_getMaxcells(maxcells)
  implicit none
  integer, intent(OUT) :: maxCells
  maxcells=max(NXB,NYB,NZB)
  !! in future if we allow use of AMReX in patch mode
  !! some implementation such as below will have to be used
  !! if(not.octree_mode) maxcells=max(largest_nxb,largest_nyb,largest_nzb)
end subroutine Grid_getMaxcells
