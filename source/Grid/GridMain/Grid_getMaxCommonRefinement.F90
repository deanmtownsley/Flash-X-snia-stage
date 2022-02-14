!!****if* source/Grid/GridMain/Grid_getMaxCommonRefinement
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
!!  Grid_getMaxCommonRefinement
!!
!! SYNOPSIS
!!
!!  call Grid_getMaxCommonRefinement(integer(IN) :: inputComm, &
!!                              integer(OUT) :: maxRefinement)
!!
!! DESCRIPTION
!!
!!  This is a simple routine to find the maximum common block refinement 
!!  level in the grid.  We may have an AMR grid in which one portion of 
!!  the grid is highly refined.  Here, it may be useful to determine the 
!!  highest block refinement level such that blocks of this level 
!!  completely cover the computational domain.
!!
!! ARGUMENTS
!!
!!  inputComm - Input MPI communicator.
!!  maxRefinement - Max common refinement level of blocks in the 
!!                  inputComm communicator.
!!
!! NOTES
!!
!!  DEV: Not yet implemented for Amrex, will abort if called.
!!
!!  Communicator argument allows us to compare a subset of processes.
!!  It also makes it explicit to the user that this routine must be 
!!  called by all processes in the passed communicator ...otherwise
!!  deadlock.
!!
!!  This routine differs from Grid_getMaxRefinement is several ways:
!!   1. Grid_getMaxCommonRefinement looks for existing LEAF blocks with the
!!      smallest refinement level (actual), while
!!      Grid_getMaxRefinement looks for blocks with the highest
!!      refinement level (either actual or potential).
!!   2. Grid_getMaxRefinement has additional optional arguments to select
!!      modes and task subsets.
!!***

subroutine Grid_getMaxCommonRefinement(inputComm, maxCommonRefinement)

#include "Simulation.h"
#include "constants.h"

  use Driver_interface, ONLY : Driver_abort
#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : lnblocks, nodetype, lrefine
#endif
#include "Flashx_mpi_implicitNone.fh"
  
  integer, intent(IN) :: inputComm
  integer, intent(OUT) :: maxCommonRefinement
  integer :: globalMaxSingleLevel, localMaxSingleLevel, lb, ierr


#if defined(FLASH_GRID_PARAMESH)
  !Find the minimum refinement level of all LEAF blocks locally, then
  !MPI_Allreduce to get global minimum.
  localMaxSingleLevel = huge(1)
  do lb = 1, lnblocks
     if (nodetype(lb) == LEAF) then
        localMaxSingleLevel = min(lrefine(lb), localMaxSingleLevel)
     end if
  end do

  call MPI_Allreduce(localMaxSingleLevel, globalMaxSingleLevel, 1, MPI_INTEGER, &
       MPI_MIN, inputComm, ierr)
  maxCommonRefinement = globalMaxSingleLevel


#elif defined(FLASH_GRID_UG)
  maxCommonRefinement = 1

#else
  maxCommonRefinement = 1
  call Driver_abort("Not yet coded")

#endif


end subroutine Grid_getMaxCommonRefinement
