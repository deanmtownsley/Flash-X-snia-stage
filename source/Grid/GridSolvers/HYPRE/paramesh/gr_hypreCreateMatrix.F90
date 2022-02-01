!!****if* source/Grid/GridSolvers/HYPRE/paramesh/gr_hypreCreateMatrix
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
!!  NAME 
!!
!!   gr_hypreCreateMatrix
!!
!!  SYNOPSIS
!!
!!   call gr_hypreCreateMatrix(integer(IN)           :: iVar,
!!                             integer(IN)           :: iFactorB,
!!                             integer(IN)           :: iFactorA,
!!                             integer(IN)           :: bcTypes(6),
!!                             real(IN)              :: bcValues(2,6),
!!                             real(IN)              :: dt,
!!                             real(IN)              :: alpha,
!!                             integer(IN)           :: blockCount,
!!                             integer(IN)           :: blockList(blockCount),
!!                             logical(IN)           :: JacobiMatrix)
!!
!!
!!  DESCRIPTION 
!!      This routine computes one of the matrices A and B, depending on the
!!      logical input argument 'JacobiMatrix':
!!          Ax = b, where A is the matrix to be inverted
!!          B = MX, where M is a matrix whose product with iVar produces RHS B.
!!
!!      A*(df/dt) + C*f = div(B*grad(f)) + D
!!      f -> Variable to be diffused.
!!      C,D are optional factors (not implemented here, the caller should add them later.)
!!
!!
!! ARGUMENTS
!! 
!!   iVar         : Variable on which the diffusion operation is performed (e.g., TEMP_VAR)
!!   iFactorB     : a factors in the equation with spatial variation.
!!   iFactorA     : a factors in the equation with spatial variation.
!!   bcTypes      : Boundary condition types.  Should be chosen from the constants
!!                  GRID_PDE_BND_* defined in Grid_interface.F90, or the special value VACUUM
!!                  defined in constants.h.
!!                  Presently OUTFLOW and VACUUM are supported, DIRICHLET less well tested.
!!   bcValues     : Values of iVar,iFactorB (!DEV: ??) on boundary (currently used for DIRICHLET and GIVENGRAD BCs).                        
!!   dt           : The time step.
!!   alpha        : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank-Nicolson
!!   blockCount   : The number of blocks in the list.   
!!   blockList    : The list of blocks on which the solution must be updated.        
!!   JacobiMatrix : TRUE computes A; FALSE computes M.
!!
!! SIDE EFFECTS
!!
!!   On return, the elements of the HYPRE matrix A that represent the second-derivative
!!   operator (div B grad) term have been defined, and it is ready for use.
!!   The gr_hypreData module variable gr_hypreMatA holds the
!!   handle for the HYPRE Solver object.
!!
!! NOTES
!!
!!   This routine does not actually 'create' the matrix object in the sense of HYPRE.
!!   It expects a matrix object already created and initialized; that is done,
!!   together with initialization of the grid object, in gr_hypreSetupGrid.
!!
!!   Currently, gr_hypreCreateMatrixFcB is called from Grid_advanceDiffusion with
!!   JacobiMatrix==.FALSE. only when the implicitness parameter theta passed to
!!   Grid_advanceDiffusion is 0. (KW 2012-12-05, corrected 2014-12-05)
!!
!! SEE ALSO
!!
!!  Grid_interface
!!***

!!REORDER(4): solnVec

subroutine gr_hypreCreateMatrix(blockCount, blockType, bcTypes, bcValues)
  
  use gr_hypreLocalInterface, ONLY: gr_hypreApplyBcToFace
  use gr_hypreData,     ONLY : gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, &
                               gr_hypreAnisoDiffusion, gr_hypreMatrixIsSetup

  use Timers_interface, ONLY : Timers_start, Timers_stop 
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
                               GRID_PDE_BND_NEUMANN,   &
                               GRID_PDE_BND_DIRICHLET

  use Grid_interface,   ONLY : Grid_getTileIterator, Grid_releaseTileIterator, Grid_getCellFaceAreas

  use Grid_tile, ONLY : Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t
 
  use Driver_interface, ONLY : Driver_abort
  
  implicit none
#include "Simulation.h"  
#include "constants.h"
#include "HYPREf.h"    
  
  !!-----------------------------------------------------------------------
  !!         ARGUMENTS
  !!-----------------------------------------------------------------------
  integer, intent(IN) :: blockCount
  integer, intent(IN) :: blockType
  integer, intent(IN) :: bcTypes(6)
  real,    intent(IN) :: bcValues(2,6)
    
  !!-----------------------------------------------------------------------
  !!         LOCAL VARIABLES.
  !!-----------------------------------------------------------------------  
  integer :: ierr
  real, dimension(MDIM)     :: del
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  integer :: datasize(MDIM), datasizeGC(MDIM)
  integer, parameter ::  mypart = 0  !! HYPRE part
  integer ::  var
  integer ::  nentries, stencil_indices(19)
  integer :: i, j, k,  lb
  integer, dimension(2,MDIM):: faces 
  real :: condimh, condiph
  real :: condjmh, condjph
  real :: condkmh, condkph
  integer :: dir, ii 
  real, allocatable :: BoxVal(:)
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor

  !! TODO : Move this out to Grid_solvePoisson level - Akash
  if (gr_hypreMatrixIsSetUp) return
 
  call Timers_start("gr_hypreCreateMatrix") 

  gr_hypreMatrixIsSetup = .TRUE.

  nullify(solnVec)
   
  nentries = NDIM*2 + 1

  do i = 1, nentries
     stencil_indices(i) = i-1
  enddo
  
  var    = 0  !! var iterator.

  lb = 1 
  call Grid_getTileIterator(itor, nodetype=blockType)
  do while(itor%isValid())

     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnVec, CENTER)

     blkLimits = tileDesc%limits
     blkLimitsGC = tileDesc%blkLimitsGC

     call tileDesc%deltas(del)
     call tileDesc%faceBCs(faces)

     datasize(1:MDIM)   = blkLimits(HIGH,1:MDIM)   - blkLimits(LOW,1:MDIM)  +1 
     datasizeGC(1:MDIM) = blkLimitsGC(HIGH,1:MDIM) - blkLimitsGC(LOW,1:MDIM)+1
     
     allocate(BoxVal(nentries*product(datasize(1:NDIM))))
    
     ii = 1
     BoxVal = 0.0

     !! We construct the matrix A (of AX=b) here, without taking into account
     !! those values in GC regions. In other words, BoxVal at
     !! ii = blkLimits(LOW,IAXIS) or blkLimits(HIGH, IAXIS) does not
     !! include condimh or condiph in this do-loop.
     !! We include these in the below by calling gr_hypreApplyBcToFace separately.
     do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
              !print*,JacobiMatrix;stop
              condimh = 0.
              condiph = 0.
              
              !! i-1,j,k
              if ((i /= blkLimits(LOW, IAXIS)) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) then                         
                 condimh = 1.
                 BoxVal(ii+1) = -condimh/(del(IAXIS)*del(IAXIS))

              end if
              
              !! i+1,j,k
              if (i /= blkLimits(HIGH, IAXIS) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) then
                 condiph = 1.
                 BoxVal(ii+2) =  -condiph/(del(IAXIS)*del(IAXIS))
              end if
              
              BoxVal(ii) =  BoxVal(ii) +  condimh/(del(IAXIS)*del(IAXIS)) + &
                                          condiph/(del(IAXIS)*del(IAXIS)) 
#if NDIM >= 2
              condjmh = 0.
              condjph = 0.

              !! i,j-1,k
              if ((j /= blkLimits(LOW, JAXIS)) .or. (faces(1,JAXIS) == NOT_BOUNDARY)) then                                    
                 condjmh = 1.
                 BoxVal(ii+3) = -condjmh/(del(JAXIS)*del(JAXIS))               
              end if
              
              !! i,j+1,k
              if ((j /= blkLimits(HIGH, JAXIS)) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) then                      
                 condjph = 1.
                 BoxVal(ii+4) = -condjph/(del(JAXIS)*del(JAXIS))                                        
              end if
              
              BoxVal(ii) =  BoxVal(ii) +  condjmh/(del(JAXIS)*del(JAXIS)) + &
                                          condjph/(del(JAXIS)*del(JAXIS))              
              
#if NDIM == 3
              condkmh = 0.
              condkph = 0.
              
              !! i,j,k-1
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then
                 condkmh = 1.                                    
                 BoxVal(ii+5) = -condkmh/(del(KAXIS)*del(KAXIS))                
              end if
              
              !! i,j,k+1
              if ((k /= blkLimits(HIGH, KAXIS)) .or. (faces(2,KAXIS) == NOT_BOUNDARY)) then                      
                 condkph = 1.
                 BoxVal(ii+6) = -condkph/(del(KAXIS)*del(KAXIS))                                        
              end if
              
              BoxVal(ii) =  BoxVal(ii) +  condkmh/(del(KAXIS)*del(KAXIS)) + &
                                          condkph/(del(KAXIS)*del(KAXIS))              
#endif             
              
#endif               
              ii = ii + nentries
              
           end do
        end do
     end do
      
     call HYPRE_SStructMatrixSetBoxValues(gr_hypreMatA, mypart, gr_hypreLower(lb,1:NDIM), & 
          gr_hypreUpper(lb,1:NDIM), var, nentries, stencil_indices(1:nentries), BoxVal(:), ierr)
 
     deallocate (BoxVal)
     
     !! We fix BoxVal at the block boundaries.
     dir = ILO_FACE
     do i = IAXIS, NDIM
        do j = LOW, HIGH
           if (faces(j,i) /= NOT_BOUNDARY) then               
              call gr_hypreApplyBcToFace(blkLimits,blkLimitsGC,mypart,var,bcTypes(dir),dir, &
                   bcValues(:,dir), del(i), gr_hypreLower(lb,:), tileDesc)
           end if
           dir = dir + 1
        end do
     end do
     
     call tileDesc%releaseDataPtr(solnVec, CENTER)
    
     lb = lb + 1 
     call itor%next()
  end do !! block
  call Grid_releaseTileIterator(itor)  
 
  !!-----------------------------------------------------------------------
  !!         THIS IS A GLOBAL CALL.
  !!-----------------------------------------------------------------------
  call HYPRE_SStructMatrixAssemble(gr_hypreMatA, ierr)    

  call Timers_stop("gr_hypreCreateMatrix") 
  
  return
  
end subroutine gr_hypreCreateMatrix
