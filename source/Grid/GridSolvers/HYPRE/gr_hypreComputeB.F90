!!****if* source/Grid/GridSolvers/HYPRE/gr_hypreComputeB
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
!!  gr_hypreComputeB
!!
!!  SYNOPSIS
!!
!!  call gr_hypreComputeB (integer, intent(IN) :: iVar
!!                         integer, intent(IN) :: iFactorB
!!                         integer, intent(IN) :: iFactorA
!!                         integer, OPTIONAL, intent(IN) :: iFactorD
!!                         real, intent(IN) :: dt
!!                         real, intent(IN) :: theta
!!                         integer, intent(IN) :: blockCount
!!                         integer,dimension(blockCount),intent(IN) :: blockList
!!                         integer, intent(IN) :: bcTypes(6)
!!                         real,    intent(IN) :: bcValues(2,6))
!!
!!
!!  DESCRIPTION 
!!   Computes the RHS of AX=B using a precomputed matrix M (stored in diff_A) such that B=MX.
!!   A MatVec product is performed to compute B and additional source terms are added 
!!   if required. A lot of the arguments are not required by this routine but are never the less
!!   provided so that B can be computed outside of HYPRE (if needed).
!!  
!!
!! ARGUMENTS
!!   iVar          : Variable on which the diffusion operatorion is performed (e.g TEMP_VAR)
!!   blockCount    : The number of blocks in the list.   
!!   blockList     : The list of blocks on which the solution must be updated.   
!!   iFactorA      :| Are factors in the equation with spatial variation. Factor D is  optional 
!!   iFactorB      :| and is generally used to represent emission in MGD. 
!!   iFactorD      :| 
!!   theta         : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   dt            : The time step (not used).
!!   bcTypes       : Presently OUTFLOW, VACUUM is supported, DIRICHLET is untested (not used).
!!   bcValues      : Values of iVar,iFactorB on boundary (DIRICHLET), (not used).
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!
!!   Stub implementation.
!!  
!!
!!***


!!REORDER(4): solnVec

subroutine gr_hypreComputeB (iVar, blockCount, blockType, bcTypes, bcValues)
 
  use Driver_interface, ONLY : Driver_abort

  use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getTileIterator, &
                               Grid_releaseTileIterator, Grid_getNumBlksFromType, &
                               Grid_getCellVolumes

  use gr_hypreLocalInterface, ONLY: gr_hypreApplyBcToFace_VecB

  use gr_hypreData,   ONLY   : gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, gr_hypreVecB, gr_hypreRefineMIN, &
                               gr_hypreUseFloor
  
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
       GRID_PDE_BND_NEUMANN,   &
       GRID_PDE_BND_DIRICHLET

  use Grid_tile, ONLY : Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t
  
  implicit none
  
#include "Simulation.h"
#include "HYPREf.h"  
#include "constants.h"
 
  integer, intent(IN) :: iVar
  integer, intent(IN) :: blockCount
  integer, intent(IN) :: blockType
  integer, intent(IN) :: bcTypes(6)
  real,    intent(IN) :: bcValues(2,6)
 
  integer :: mylevel, mypart, var, ii,i,j,k, ierr, dir
  real, allocatable :: RHSVal(:)
  integer :: datasize(MDIM)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer :: lb
  real, allocatable :: cellVolumes(:,:,:)
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  logical :: mask(NUNK_VARS), savedUseFloor

  real, dimension(MDIM) :: del
  integer, dimension(2,MDIM):: faces

  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor

  nullify(solnVec)

  var = 0
  mypart = 0

  lb = 1
  call Grid_getTileIterator(itor, nodetype=blockType)
  do while(itor%isValid())

     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnVec, CENTER)

     blkLimits   = tileDesc%limits
     blkLimitsGC = tileDesc%blkLimitsGC

     call tileDesc%deltas(del)
     call tileDesc%faceBCs(faces)
      
     mylevel = tileDesc%level

     mypart = mylevel - gr_hypreRefineMIN
     
     datasize(1:MDIM) = blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM)+1         
         
     allocate(RHSVal(product(dataSize(1:NDIM))))     
 
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
              
              ii = (k - blkLimits(LOW,KAXIS)  + 1)                             +  &
                   (j - blkLimits(LOW,JAXIS))*dataSize(KAXIS)                  +  &
                   (i - blkLimits(LOW,IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)  
              
              RHSVal(ii) = -solnVec(iVar,i,j,k)
              
           end do
        end do
     end do
              
     call HYPRE_SStructVectorSetBoxValues(gr_hypreVecB, mypart, gr_hypreLower(lb, 1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, RHSVal(:), ierr)  
     
 
     deallocate (RHSVal)

     !! We fix RHSVal at the block boundaries.
     dir = ILO_FACE
     do i = IAXIS, NDIM
        do j = LOW, HIGH
           if (faces(j,i) /= NOT_BOUNDARY) then               
              call gr_hypreApplyBcToFace_VecB(blkLimits,blkLimitsGC,mypart,var,bcTypes(dir),dir, &
                   bcValues(:,dir), del(i), gr_hypreLower(lb,:), tileDesc)
           end if
           dir = dir + 1
        end do
     end do
 
     lb = lb + 1

     call tileDesc%releaseDataPtr(solnVec, CENTER)
     call itor%next()

  end do
  call Grid_releaseTileIterator(itor)  

  call HYPRE_SStructVectorAssemble(gr_hypreVecB, ierr)  
 
end subroutine gr_hypreComputeB
