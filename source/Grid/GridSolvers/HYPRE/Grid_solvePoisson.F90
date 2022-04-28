!!****if* source/Grid/GridSolvers/HYPRE/Grid_solvePoisson
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
!!  Grid_solvePoisson
!!
!! SYNOPSIS
!!
!!  call Grid_solvePoisson(integer(IN) :: iSoln,
!!                         integer(IN) :: iSrc, 
!!                         integer(IN) :: bcTypes(6),
!!                         real(IN)    :: bcValues(2,6),
!!                         real(INOUT) :: poisfact)
!!
!! DESCRIPTION
!!
!!   Driver routine for poisson solvers in the Grid
!!
!!
!! ARGUMENTS
!!
!!  iSoln    - the index for the solution variable (potential when used for self-gravity)
!!  iSrc     - the index of the source variable (density when used for self-gravity)
!!  bcTypes  - the boundary condition type; only the first entry is used.
!!             Only the first 2*NDIM elements are significant. They are interpreted
!!             in the order (X left, X right, Y left, Y right, Z left, Z right).
!!             Valid values are:
!!               GRID_PDE_BND_PERIODIC (1)
!!               GRID_PDE_BND_DIRICHLET (2) (homogeneous or constant Dirichlet)
!!               GRID_PDE_BND_NEUMANN (3) (homogeneous or constant Neumann)
!!               GRID_PDE_BND_ISOLATED (0)
!!           !DEV: requested bcTypes ignored, always uses GRID_PDE_BND_NEUMANN
!!  bcValues - the values to boundary conditions, currently not used (treated as 0)
!!  poisfact      - scaling factor to be used in calculation
!!
!! NOTES
!!
!!  The symbols listed above for bcTypes are declared as FORTRAN PARAMETERS in
!!  the module Grid_interfaces.  Code using this interface should refer to that
!!  module with a USE statement, like this:
!!
!!    use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
!!       GRID_PDE_BND_ISOLATED, GRID_PDE_BND_DIRICHLET, &
!!       Grid_solvePoisson
!!  
!!  Most implementations only support a limited subset of boundary condition types.
!!  Some implementations require that all significant elements of bcTypes are the same.
!!  (That is always the case for GRID_PDE_BND_ISOLATED.)
!!  Even if an implementation supports combinations of different boundary conditions
!!  on different sides of the domain, the types at left and right sides for the same
!!  axis direction will usually have to be the samme.
!!
!!  Support in some implementations provided with FLASH4:
!!
!!   GridSolvers/Multipole:                  GRID_PDE_BND_ISOLATED                   
!!   GridSolvers/Multigrid (simple):         GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.),
!!                                            GRID_PDE_BND_NEUMANN(hom.)(?)
!!                                            (same type in all directions),
!!                                            or GRID_PDE_BND_ISOLATED
!!                                           (requires Paramesh as Grid with NBlockX==NBlockY==NBlockZ==1)
!!   GridSolvers/Pfft:                       GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET(hom.), 
!!                                            GRID_PDE_BND_NEUMANN(hom.)
!!                                            in various combinations,
!!                                            depending on GridSolvers/Pfft subdirectory 
!!                                            (i.e. implementation configured in)
!!                                           (requires UG in pencil shape or Paramesh as Grid)
!!   GridSolvers/Multigrid hybrid with Pfft: GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET, 
!!                                            GRID_PDE_BND_NEUMANN in various combinations,
!!                                            or GRID_PDE_BND_ISOLATED
!!                                           (requires Paramesh as Grid)
!!   GridSolvers/HYPRE:                      currently GRID_PDE_BND_NEUMANN only
!!   
!!***

subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact)

  use Grid_data,        ONLY : gr_meshMe, gr_meshcomm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abort
  use gr_hypreLocalinterface,  ONLY : gr_hypreCreateMatrix, gr_hypreComputeB,    &
                               gr_hypreGridStatus, gr_hypreSetIniGuess, &
                               gr_hypreUpdateSoln, gr_hypreSolve
  use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getTileIterator, &
                               Grid_releaseTileIterator, Grid_getNumBlksFromType, &
                               Grid_getCellVolumes

  use gr_hypreData,   ONLY   : gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, gr_hypreVecB, gr_hypreRefineMIN, &
                               gr_hypreUseFloor, gr_hypreMatrixIsSetup
 
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
       GRID_PDE_BND_NEUMANN,   &
       GRID_PDE_BND_DIRICHLET

  use Grid_tile, ONLY : Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t
 
  implicit none

#include "Simulation.h"
#include "HYPREf.h"  
#include "constants.h"
  
  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes(6)
  real, intent(in)       :: bcValues(2,6)
  real, intent(inout)    :: poisfact !DEV: NOT intent(IN) because some implementation actually changes it? - KW  
  
  logical :: mask(NUNK_VARS), savedUseFloor

  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor

  integer :: blockCount, blockType

  integer TA(2),count_rate
  real :: ET
 
!!$    character(len=32) :: matfile
  
  call Timers_start("Grid_solvePoisson")    


  blockType  = LEAF
  call Grid_getNumBlksFromType(blockType,blockCount)

  !!-----------------------------------------------------------------------
  !!     1.  Do we need to reset HYPRE grid ?, has the underlying AMR
  !!         mesh been modified ?, is this the first call to HYPRE ?
  !!         Only in AMR is this routine meaningful.
  !!-----------------------------------------------------------------------
  call gr_hypreGridStatus (blockCount,blockType)

  !!-----------------------------------------------------------------------
  !!     1. Compute RHS
  !!
  !!----------------------------------------------------------------------- 
  call gr_hypreComputeB (iSrc, blockCount, blockType, bcTypes, bcValues) 

  !!-----------------------------------------------------------------------
  !!     2. Set initial guess for solver typically 
  !!        iVar at previous time step is used.
  !!-----------------------------------------------------------------------
  call gr_hypreSetIniGuess (iSoln,blockCount,blockType)      

  !!-----------------------------------------------------------------------
  !!     3. Create A Matrix
  !!-----------------------------------------------------------------------
  if (gr_meshMe .eq. 0 .and. gr_hypreMatrixIsSetup) then
     write(*,*) 'Using HYPRE matrix from previous time step'
  else if (gr_meshMe .eq. 0.0) then
     write(*,*) 'Reconstructing HYPRE matrix'
  end if 
  call gr_hypreCreateMatrix (blockCount, blockType, bcTypes, bcValues)
 
  !!-----------------------------------------------------------------------
  !!     4. Solve AX = B
  !!-----------------------------------------------------------------------
  if (gr_meshMe .eq. 0) CALL SYSTEM_CLOCK(TA(1),count_rate)  
  call gr_hypreSolve ()
  if (gr_meshMe .eq. 0) then
     CALL SYSTEM_CLOCK(TA(2),count_rate)
     ET=REAL(TA(2)-TA(1))/count_rate
     write(*,*) 'HYPRE Poisson Solver time = ',ET,' sec.'
  endif

  !!-----------------------------------------------------------------------
  !!     5. Update unk variable using X.
  !!-----------------------------------------------------------------------
  savedUseFloor = gr_hypreUseFloor
  gr_hypreUseFloor = .FALSE.    ! Do not floor solution of Poisson problem.
  call gr_hypreUpdateSoln (iSoln, blockCount, blockType)  
  gr_hypreUseFloor = savedUseFloor    ! rstore previous value of flag
  
  call Timers_stop("Grid_solvePoisson") 
  
  return
end subroutine Grid_solvePoisson
