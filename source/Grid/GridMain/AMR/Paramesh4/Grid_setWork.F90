!!****f* source/Grid/GridMain/AMR/Paramesh4/bittree/Grid_setWork
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
!!  Grid_setWork
!!
!! SYNOPSIS
!!
!!  Grid_setWork(Grid_tile_t(IN)       :: tileDesc
!!               real(IN)              :: work   )
!!               integer(IN),optional  :: mode   )
!!  
!! DESCRIPTION 
!!  Sets the work value for a given block. Can pass an operation name (defined in
!!  constants.h) to apply an operation to the work instead of overwriting the
!!  previous value. Used for uneven load distributions when using Paramesh.
!!
!!  The optional mode argument controls the behavior of updating work. The user
!!  can pass values GRIDOP_{SET,ADD,SUB,MLT,DIV,AVG,MAX,MIN} which replace the old value
!!  add to the old value, etc. These constants are defined in constants.h
!!
!! ARGUMENTS
!!  tileDesc - block tile.
!!  work - desired weight in the work distribution.
!!  mode - (optional) If supplied, specifies the operation to apply the given work.
!!         Default = GRIDOP_SET (replace old value with passed value).
!!
!!***
#include "constants.h"

      subroutine Grid_setWork(tileDesc, work, mode)
      use Grid_tile, ONLY : Grid_tile_t
      use tree, only : work_block, nodetype, &
                       gr_btCustomWork, gr_btWorkBoundsPar, &
                       gr_btWorkBoundsLeaf
      use Driver_interface, only : Driver_abort

      implicit none
      type(Grid_tile_t),intent(in) :: tileDesc
      real,intent(in)              :: work
      integer,intent(in),optional  :: mode

      integer           :: opmode, blkid
      real              :: work_new, block_default
      integer,parameter :: r8 = kind(1.0)
      real,parameter    :: eps=1.0e-100_r8

#ifdef FLASH_DEBUG_AMR
      if(.NOT.gr_btCustomWork) &
        call Driver_abort( &
                "Grid_setWork: Trying to set work array, &
                &but simulation is not configured to sort via custom &
                &work values at regridding. Use `gr_btCustomWork &
                &= True` in your par file for desired results." )
#endif

!-----Get block ID
      blkid = tileDesc%id

!-----Set operation mode (default = replace previous value)
      opmode = GRIDOP_SET
      if(present(mode)) opmode = mode

!-----Determine new value for work
      select case (opmode)
        case (GRIDOP_SET)
          work_new = work
        case (GRIDOP_ADD)
          work_new = work_block(blkid) + work
        case (GRIDOP_SUB)
          work_new = work_block(blkid) - work
        case (GRIDOP_MLT)
          work_new = work_block(blkid) * work
        case (GRIDOP_DIV)
          if (abs(work).lt.eps) then
            call Driver_abort("Grid_setWork: &
                    &tried to divide current work by zero")
          else
            work_new = work_block(blkid) / work
          end if
        case (GRIDOP_AVG)
          work_new = (work_block(blkid) + work) / 2.0
        case (GRIDOP_MAX)
          work_new = max(work_block(blkid),work)
        case (GRIDOP_MIN)
          work_new = min(work_block(blkid),work)
        case default
          call Driver_abort("Unknown mode passed to Grid_setWork.")
      end select

!-----Enforce bounds on work
      if(nodetype(blkid).eq.1) then
        work_new = max(work_new, gr_btWorkBoundsLeaf(LOW))
        work_new = min(work_new, gr_btWorkBoundsLeaf(HIGH))
      else
        work_new = max(work_new, gr_btWorkBoundsPar(LOW))
        work_new = min(work_new, gr_btWorkBoundsPar(HIGH))
      end if

!-----Store new value in work_block
      work_block(blkid) = work_new

      return
      end subroutine Grid_setWork
