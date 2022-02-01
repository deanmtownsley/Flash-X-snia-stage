!!****if* source/Grid/GridMain/AMR/Amrex/Grid_getLocalNumBlks
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
!!  Grid_getLocalNumBlks
!!
!! SYNOPSIS
!!
!!  call Grid_getLocalNumBlks(integer(OUT) :: numBlocks)
!!
!! DESCRIPTION
!!  Get the number of local blocks on a processor
!!
!! ARGUMENTS
!!  numBlocks : The number of blocks currently in use on myProcessor
!!
!! NOTES
!!  There should be a better (more efficient) way to do this!
!!  The current version iterates through a loop and counts interations.
!!  There should be (and maybe already is) a better way to get the
!!  information directly from AMReX.
!!
!!  The blocks counted include LEAF as well as covered blocks.
!!  This is consistent with the functionality of the PARAMESH version, and used
!!  in this way in the friendly IO unit.
!!
!!  An alternative version Grid_getLocalNumLeafBlks is coded below,
!!  but currently (2022-01-25) not yet made public via Grid_interface.
!!***

#include "constants.h"

subroutine Grid_getLocalNumBlks(numBlocks)
  use Grid_interface, ONLY : Grid_getTileIterator, &
                             Grid_releaseTileIterator
  use Grid_iterator,  ONLY : Grid_iterator_t

  implicit none

  integer, intent(OUT) :: numBlocks

  type(Grid_iterator_t) :: itor

  numBlocks = 0

  call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.)
  do while (itor%isValid())
     numBlocks = numBlocks + 1
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
end subroutine Grid_getLocalNumBlks

subroutine Grid_getLocalNumLeafBlks(numBlocks)
  use Grid_interface, ONLY : Grid_getTileIterator, &
                             Grid_releaseTileIterator
  use Grid_iterator,  ONLY : Grid_iterator_t

  implicit none

  integer, intent(OUT) :: numBlocks

  type(Grid_iterator_t) :: itor

  numBlocks = 0

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while (itor%isValid())
     numBlocks = numBlocks + 1
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
end subroutine Grid_getLocalNumLeafBlks

