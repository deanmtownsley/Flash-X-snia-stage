!!****if* source/Grid/GridMain/paramesh/gr_xyzToBlock
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
!!  gr_xyzToBlock
!!
!! SYNOPSIS
!!
!!  call gr_xyzToBlock( real (IN)  :: xyz(MDIM),
!!                    integer(OUT) :: procID,
!!                    integer(OUT) :: blkID)
!!  
!! DESCRIPTION 
!!  
!!  This routine returns the identity of the block, on
!!  which the specified physical location falls. It also
!!  returns the ID of the processor on which the block is
!!  residing.
!!
!!  If the physical location lies outside the domain
!!  boundaries, then procID=NONEXISTENT will be returned;
!!  boundary periodicity is not considered.
!!
!! ARGUMENTS 
!!
!! NOTES
!!
!!  This implementation requires BITTREE.
!!***

#include "constants.h"
#include "Simulation.h"

Subroutine gr_xyzToBlock(xyz, procID, blkID)
  use gr_interface, ONLY : gr_xyzToBlockLevel
  use Grid_data, ONLY : gr_meshNumProcs
  use Driver_interface, ONLY : Driver_abort
  use tree, ONLY : lrefine_max
#ifdef BITTREE
  use bittree, only : gr_btIdentify
#endif

  implicit none

  real, dimension(MDIM),intent(IN) :: xyz
  integer, intent(OUT) :: procID
  integer, intent(OUT) :: blkID

#ifdef BITTREE
  integer, dimension(MDIM) :: ijk
#endif
  integer :: lev, proc, blk

  lev = lrefine_max

#ifdef BITTREE  
  call gr_xyzToBlockLevel(lev, xyz, ijk)
  call gr_btIdentify(gr_meshNumProcs, lev, ijk, proc, blk)
  if (lev < 1) then       !if the point was outside of the domain...
     proc = NONEXISTENT         !make sure we communicate nonexistence to the caller.
  end if
#else
  proc = -size(xyz); blk= -1  !avoid compiler warnings about uninitialized variables below
  call Driver_abort("gr_xyzToBlock works only when bittree is enabled")
#endif

  procID=proc
  blkID=blk
End Subroutine gr_xyzToBlock
