!!****if* source/Grid/GridMain/AMR/Grid_getTileIterator
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
!!  Grid_getTileIterator
!!
!! SYNOPSIS
!!  Grid_getTileIterator(Grid_iterator_t(OUT)  :: itor,
!!                       integer(IN)           :: nodetype,
!!                       integer(IN), optional :: level,
!!                       logical(IN), optional :: tiling,
!!                       integer(IN), optional :: nthreads)
!!  
!! DESCRIPTION 
!!  Construct an iterator for walking across a specific subset of blocks or
!!  tiles within the current octree structure.  The iterator is already
!!  set to the first matching leaf block/tile.
!!
!!  Once finished, the iterator should be destroyed with
!!  Grid_releaseTileIterator
!!
!! ARGUMENTS 
!!  itor     - the requested block/tile iterator
!!  nodetype - the class of blocks to iterate over.  The options for AMReX are
!!                - ALL_BLKS - all blocks
!!                - LEAF     - only leaf blocks
!!  level    - iterate only over leaf blocks/tiles located at this level of
!!             refinement.  A level value of UNSPEC_LEVEL is equivalent to
!!             omitting this optional argument.
!!  tiling   - an optional optimization hint.  If TRUE, then the iterator will
!!             walk across all associated blocks on a tile-by-tile basis *if*
!!             the implementation supports this feature.  If a value is not
!!             given, is FALSE, or the implementation does not support tiling,
!!             the iterator will iterate on a block-by-block basis.
!!  nthreads - an optional argument to change the behavior of the iterator
!!             in an OpenMP parallel region from the default.
!!             By default, if Grid_getTileIterator is called in an active
!!             parallel region, the iteration space will be divided up so
!!             that each thread will end up executing an approximately equal
!!             number of iterations, assuming that all threads will call
!!             the iterator's next() method repeatedly until isValid() returns
!!             FALSE. This is equivalent to what the MFIter of AMReX does
!!             when it is not in 'dynamic' mode, and is also very similar
!!             to the static scheduling of do / for looks by OpenMP.
!!             The described default behavior in an active parallel region
!!             is achieved by effectively modifying the default serial
!!             behavior in two ways:
!!             o Each thread begins iterating at a thread-specific starting
!!               index;
!!             o when the next() method is called in a thread, the index
!!               position advances not by 1 but by a number n.
!!             This number n is by default set to the return value of
!!             omp_get_num_threads(). If the optional argument 'nthreads'
!!             is present, n is set to this value. Moreover, if 'nthreads'
!!             is present, then the thread-specific starting index mentioned
!!             above will be set to 0 (zero-based convention) if it would
!!             otherwise exceed nthreads-1.
!!             This argument will only have and effect if all of the following
!!             are true:
!!             o The code is built with OpenMP support (-fopenmp copmiuler flag
!!               or similar, depending of compiler used). (This implies that
!!               preprocessor symbol _OPENMP is defined.)
!!             o Grid_getTileIterator is called in an active parallel region
!!               (i.e., omp_get_num_threads() return a value > 1).
!!
!! NOTES
!!
!!   The described default behavior in an active parallel region
!!   makes sense in the following scenario:
!!             o The iterator object 'itor' is private to threads;
!!             o all threads active at the time of iterator construction
!!               will participate in iterations.
!!   There is, however, an alternative scenario where those
!!   characteristics do not hold true, and instead:
!!             o The iterator object 'itor' is shared by threads
!!               (or perhaps only the version belonging to one designated
!!               thread matters);
!!             o only one thread (or maybe a subset) iterates; other threads
!!               may be participating in computational work via OpenMP tasking
!!               constructs or by other means, but this division of labor,
!!               if it exists, is not visible to the iterator.
!!
!! SEE ALSO
!!  Grid_releaseTileIterator.F90
!!
!!***

#include "constants.h"

subroutine Grid_getTileIterator(itor, nodetype, level, tiling, tileSize, nthreads)
  use Grid_iterator, ONLY : Grid_iterator_t, build_iterator
  use Grid_data,     ONLY : gr_enableTiling, &
                            gr_useTiling,    &
                            gr_tileSize

  implicit none

  type(Grid_iterator_t), intent(OUT)          :: itor
  integer,               intent(IN)           :: nodetype
  integer,               intent(IN), optional :: level
  logical,               intent(IN), optional :: tiling
  integer,               intent(IN), optional :: tileSize(1:MDIM)
  integer,               intent(IN), optional :: nthreads

  integer :: myLevel
  logical :: myTiling
  integer :: myTileSize(1:MDIM)

  myLevel = UNSPEC_LEVEL
  if (present(level)) then
    myLevel = level
  end if

  myTiling = .FALSE.
  if (gr_enableTiling) then
     if (present(tiling)) then
        myTiling = tiling
     else
        myTiling = gr_useTiling
     end if
  end if

  if (present(tileSize)) then
    myTileSize = tileSize 
  else
    myTileSize = gr_tileSize
  end if

  call build_iterator(itor, nodetype, myLevel, myTiling, myTileSize, nthreads=nthreads)
end subroutine Grid_getTileIterator

