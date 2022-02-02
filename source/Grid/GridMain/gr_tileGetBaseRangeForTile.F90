#include "constants.h"
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
#include "Simulation.h"

subroutine gr_tileGetBaseRangeForTile(range, tileID)
  use Grid_getBlkIndexLimits_mod, ONLY : Grid_getBlkIndexLimits
  use gr_tilePolicyData,ONLY: nxt=>gr_tileNxt,nyt=>gr_tileNyt,nzt=>gr_tileNzt
  implicit none

#ifdef FIXEDBLOCKSIZE
  integer, parameter :: nxb=NXB,nyb=NYB,nzb=NZB
#else
  integer :: nxb,nyb,nzb
#endif
  integer, intent(out) :: range(2,MDIM)
  integer, intent(IN) :: tileID

  integer :: blockID
  integer :: blkLimitsGcUnused(2,MDIM)
  integer :: tx,ty,tz

  if (tileID < 100000) then
     call Grid_getBlkIndexLimits(tileID, range, blkLimitsGcUnused, CENTER)
  else
     blockID = tileID / 100000
     call Grid_getBlkIndexLimits(blockID, range, blkLimitsGcUnused, CENTER)
#ifndef FIXEDBLOCKSIZE
     nxb = range(HIGH,IAXIS)-range(LOW,IAXIS)+1
     nyb = range(HIGH,JAXIS)-range(LOW,JAXIS)+1
     nzb = range(HIGH,KAXIS)-range(LOW,KAXIS)+1
#endif
     call gr_tile2txtytz(tileID,tx,ty,tz)
     range(2,IAXIS) = min(range(2,IAXIS), range(1,IAXIS)+nxb*(tx+1)/nxt-1)
     range(1,IAXIS) =     range(1,IAXIS) + nxb*tx/nxt
#if (NDIM > 1)
     range(2,JAXIS) = min(range(2,JAXIS), range(1,JAXIS)+nyb*(ty+1)/nyt-1)
     range(1,JAXIS) =     range(1,JAXIS) + nyb*ty/nyt
#endif
#if (NDIM > 2)
     range(2,KAXIS) = min(range(2,KAXIS), range(1,KAXIS)+nzb*(tz+1)/nzt-1)
     range(1,KAXIS) =     range(1,KAXIS) + nzb*tz/nzt
#endif

  end if
end subroutine gr_tileGetBaseRangeForTile
