#include "constants.h"

!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! A Milhoja-specific implementation of this routine.  Please refer
!! to the documentation in this routine's stub for general interface
!! information.
!!
!! @todo There should be a better (more efficient) way to do this!
!!  The current version iterates through a loop and counts iterations.
!!  There should be (and maybe already is) a better way to get the
!!  information directly from the Milhoja grid backend.
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

