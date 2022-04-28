!!****if* source/Grid/GridMain/Grid_getMaxRefinement
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
!!  Grid_getMaxRefinement
!!
!! SYNOPSIS
!!
!!  call Grid_getMaxRefinement(integer(OUT) :: maxRefinement,
!!                    OPTIONAL,integer(IN)  :: mode,
!!                    OPTIONAL,integer(IN)  :: scope,
!!                    OPTIONAL,integer(IN)  :: inputComm
!!                              )
!!
!! DESCRIPTION
!!
!!  This routine returns the maximum block refinement level in the grid.  
!!  Depending on the mode used (and modified by other optional arguments),
!!  the returned value represents the maximum refinement level that
!!  is either allowed, or is currently realized anywhere on the grid
!!  (or a subset of the grid).
!!
!!  We may have an AMR grid in which one portion of the grid is highly
!!  refined.  Here, it may be useful to determine globally the highest block
!!  refinement level that occurs, or can potentially occur during the
!!  simulation.
!!
!! ARGUMENTS
!!
!!  maxRefinement - Maximum refinement level returned.
!!
!!  inputComm - Input MPI communicator, only used if mode=4 and scope=2.
!!              Default - none.
!!  mode      - 1 for lrefine_max,
!!              2 for gr_maxRefine (if used) or lrefine_max,
!!              3 for gr_maxRefine,
!!              4 for existing blocks.
!!              Default is 3.
!!  scope     - scope only used if mode=4;
!!              1 for local to this MPI task,
!!              2 for tasks in inputComm,
!!             [3 for mesh communicator,]
!!             [4 for MPI_COMM_WORLD].
!!              Default is 3.
!!
!! NOTES
!! 
!!  Communicator argument allows us to compare a subset of processes.
!!  It also makes it explicit to the user that this routine must be 
!!  called by all processes in the passed communicator ... otherwise
!!  deadlock. This applies only for modes that require communication.
!!
!!  For a uniform grid implementation, the returned level will always
!!  be 1.
!!
!!  This routine differs from Grid_getMaxCommonRefinement is several ways:
!!   1. Grid_getMaxCommonRefinement looks for existing LEAF blocks with the
!!      smallest refinement level (actual), while
!!      Grid_getMaxRefinement looks for blocks with the highest
!!      refinement level (either actual or potential).
!!   2. Grid_getMaxRefinement has additional optional arguments to select
!!      modes and task subsets.
!!
!!   For the Amrex implementation, not all values for scope are supported.
!!   In the case mode=4 it is up to the AMReX function amrex_get_finest_level()
!!   whether actual communication takes place.
!!
!!   In the Paramesh4 implementation, for the common case mode=4,scope=3
!!   a value from the module variable gr_finestExistingLevel may be
!!   returned rather than performing actual MPI communication; this
!!   value is assumed to be a valid cached value set to the finest
!!   existing level if it is > 0. However, this subroutine does not
!!   modify gr_finestExistingLevel; such modification should be done
!!   outside of this routine whenever the refinement pattern changes.
!!***

subroutine Grid_getMaxRefinement(maxRefinement, mode, scope, inputComm)

#include "Simulation.h"
#include "constants.h"

  use Driver_interface, ONLY : Driver_abort
  use Grid_data, ONLY : gr_meshComm, gr_meshMe
#ifdef FLASH_GRID_PARAMESH
  use Grid_data, ONLY : gr_maxRefine, gr_enforceMaxRefinement
  use gr_specificData, ONLY : gr_finestExistingLevel
  use tree, ONLY : lnblocks, lrefine_min,lrefine_max, lrefine
#endif
#ifdef FLASH_GRID_AMREX
  use amrex_amrcore_module, ONLY : amrex_max_level, amrex_get_finest_level
  use Grid_data, ONLY :  gr_maxRefine, gr_enforceMaxRefinement
#endif

#include "Flashx_mpi_implicitNone.fh"

  integer, intent(IN), OPTIONAL :: mode, scope
  integer, intent(IN), OPTIONAL :: inputComm
  integer, intent(OUT) :: maxRefinement

  integer :: ierr
  integer :: myMode, myScope, outLevel, localOutLevel
  integer :: comm
  logical :: needReduce

  if (present(mode)) then
     myMode = mode
  else
     myMode = 3
  end if

  if (present(scope)) then
     myScope = scope
  else if (myMode==4) then
     myScope = 3
  else
     myScope = -1
  end if

  if (present(inputComm) .AND. myScope==2) then
     comm = inputComm
  else if (myScope==2 .AND. myMode==4) then
     call Driver_abort('Grid_getMaxRefinement: MPI communicator not specified!')
  else if (myScope==3) then
     comm = gr_meshComm
  else if (myScope==4) then
     comm = MPI_COMM_WORLD
  else
     comm = -1
  end if

  needReduce = (myMode == 4 .AND. myScope > 1)

#ifdef FLASH_GRID_UG

  maxRefinement = 1

#elif defined FLASH_GRID_AMREX
  ! AMReX uses 0-based level indexing/FLASH uses 1-based
  
  ! Assume mode 1
  maxRefinement = amrex_max_level + 1

  if      ((myMode == 2) .AND. gr_enforceMaxRefinement)  then
    maxRefinement = MIN(maxRefinement, gr_maxRefine)
  else if  (myMode == 3) then
    maxRefinement = gr_maxRefine
  else if  (myMode == 4) then
    maxRefinement = amrex_get_finest_level() + 1
  end if

  if (present(scope)) then
     if (myMode == 4 .AND.  myScope .NE. 3) then
        call Driver_abort("[Grid_getMaxRefinement] scope not coded yet for AMReX")
     end if
  else if (present(inputComm)) then
    call Driver_abort("[Grid_getMaxRefinement] inputComm not coded yet for AMReX")
  end if

#else

# if defined(FLASH_GRID_PARAMESH)
  if (myMode == 4 .AND. myScope == 3 .AND. .NOT.present(inputComm) &
       .AND. gr_finestExistingLevel > 0) then
     maxRefinement = gr_finestExistingLevel
     return
  end if
# endif

  if (myMode .LE. 2) then
     outLevel = lrefine_max
  else
     outLevel = huge(1)
  end if

  if (myMode==2 .AND. gr_enforceMaxRefinement) then
     outLevel = min(outLevel, gr_maxRefine)
  end if

  if (myMode==3) then
     outLevel = gr_maxRefine
  end if

  if (myMode==4) then
     outLevel = maxval(lrefine(1:lnblocks))
  end if

  localOutLevel = max(0,min(lrefine_max,outLevel))

  if (.NOT. needReduce) then
#ifdef DEBUG_GRID
     print*,'[Grid_getMaxRefinement]: returning',localOutLevel,' on',gr_meshMe
#endif
     maxRefinement = localOutLevel
     return
  end if

! The following code should be reached only if defined(FLASH_GRID_PARAMESH)

  call MPI_Allreduce(localOutLevel, outLevel, 1, MPI_INTEGER, &
       MPI_MAX, comm, ierr)
#ifdef DEBUG_GRID
  print*,'[Grid_getMaxRefinement]: returning',outLevel,' on',gr_meshMe
#endif
  maxRefinement = outLevel

#endif


end subroutine Grid_getMaxRefinement
