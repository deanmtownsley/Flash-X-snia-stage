!!****if* source/Grid/GridMain/AMR/Amrex/Grid_getNumBlksFromType
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
!!  Grid_getNumBlksFromType
!!
!! SYNOPSIS
!!
!!  call Grid_getNumBlksFromType(integer(IN) :: blockType, integer(OUT) :: numBlocks)
!!  
!! DESCRIPTION 
!!  Get the number of local blocks on a processor based on type 
!!
!! ARGUMENTS
!!  blockType : The nodetype of blocks for which numBlocks is requested (ALL_BLKS or LEAF)
!!  numBlocks : The number of blocks of the requested type currently on myProcessor.
!!
!!***

#include "constants.h"

subroutine Grid_getNumBlksFromType(blockType,numBlocks)
  use Grid_interface, ONLY : Grid_getTileIterator, &
                             Grid_releaseTileIterator
  use Grid_iterator,  ONLY : Grid_iterator_t
  use Driver_interface, only: Driver_abort


  implicit none

  integer,intent(in)  :: blockType
  integer,intent(out) :: numBlocks

  type(Grid_iterator_t) :: itor

  if ( (blockType.ne.ALL_BLKS) .and. (blockType.ne.LEAF)) then
    call Driver_abort("Invalid nodetype passed to Grid_getNumBlksFromType. Valid options: ALL_BLKS, LEAF.")
  endif

  numBlocks = 0
  call Grid_getTileIterator(itor, blockType, tiling=.FALSE.)
  do while (itor%isValid())
    numBlocks = numBlocks + 1
    call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

end subroutine Grid_getNumBlksFromType
