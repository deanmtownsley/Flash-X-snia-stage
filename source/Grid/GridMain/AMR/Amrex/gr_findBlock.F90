!!****if* source/Grid/GridMain/Chombo/gr_findBlock
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
!!  gr_ptFindBlock
!!
!! SYNOPSIS
!!
!!  call gr_findBlock(integer(in)    :: blkList(blkCount),
!!                    integer(in)    :: blkCount,
!!                    real   (in)    :: pos(MDIM),
!!                    integer(INOUT) :: blockID)
!!
!! DESCRIPTION
!!
!!   Given a point in the domain, this routine finds if the
!!   point lies within one of the blocks on this process. If
!!   such a block is found, its ID is returned, otherwise
!!   blockID is set to NONEXISTENT
!!
!! ARGUMENTS
!!
!!   blkList  : The list of blocks to be examined
!!   blkCount : the number of blocks in the list
!!   pos      : the coordinates of the point
!!   blockID  : the identity of the block if found, NONEXISTENT otherwise
!!              NONEXISTENT is a constant defined in constants.h
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine gr_findBlock(blkList,blkCount,pos,blockID)
  use Driver_interface, ONLY : Driver_abort
  implicit none
  integer,intent(IN) :: blkCount
  integer,dimension(blkCount),intent(IN) :: blkList
  real,dimension(MDIM),intent(IN) :: pos
  integer,intent(INOUT) :: blockID

  ! DEV: TODO Implement this with blockDesc if needed
  call Driver_abort("[gr_findBlock] Not yet implemented for AMReX")
end subroutine gr_findBlock

