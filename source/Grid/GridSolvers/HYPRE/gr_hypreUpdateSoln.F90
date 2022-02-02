!!****if* source/Grid/GridSolvers/HYPRE/gr_hypreUpdateSoln
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
!!  gr_hypreUpdateSoln
!!
!!  SYNOPSIS
!!
!!  call gr_hypreUpdateSoln (integer,intent(IN) :: iVar,
!!                           integer,intent(IN) :: blockCount,
!!                           integer,intent(IN) :: blockList (blockCount))
!!
!!  DESCRIPTION 
!!      This routine updates solution after solve (diffusion operation AX=B).
!!
!! ARGUMENTS
!!   iVar       : Variable on which the diffusion operation is performed (e.g., TEMP_VAR)
!!   blockCount : The number of blocks in the list.   
!!   blockList  : The list of blocks on which the solution must be updated.  
!!
!!
!! SIDE EFFECTS
!!
!!  Also sets the HYPRE vector for B to 0.
!!  
!! NOTES
!!
!!   Uses HYPRE library.
!!   Solution is floored if gr_hypreUseFloor is true.
!!
!!***

!!REORDER(4): solnVec


subroutine gr_hypreUpdateSoln (iVar, blockCount, blockType)
  
  use gr_hypreData,     ONLY : gr_hypreVecX, gr_hypreLower, &
       gr_hypreRefineMIN, gr_hypreUpper, gr_hypreVecB, gr_hypreUseFloor, gr_hypreFloor, &
       gr_hypreSolnIsDelta
  use gr_solversData,   ONLY : dbgContextSolvers => gr_solversDbgContext
  use Timers_interface, ONLY : Timers_start, Timers_stop    
 
  use Grid_interface,   ONLY : Grid_getTileIterator, Grid_releaseTileIterator

  use Grid_tile, ONLY : Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t
  
  implicit none
  
#include "Simulation.h"  
#include "constants.h"
#include "HYPREf.h"  
  
  integer,intent(IN) :: iVar
  integer,intent(IN) :: blockCount
  integer,intent(IN) :: blockType
  
  !! LOCAL VARIABLES
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits  
  real :: values(1)
  integer :: blockID,part,level,var
  integer :: i, j, k, lb, ierr, pos(NDIM)
  integer :: ii
  real, allocatable, dimension(:) :: BoxVal
  integer :: datasize(MDIM)
  integer :: component
  logical :: skipUpdate
 
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor
      
  call Timers_start("gr_hypreUpdateSoln")    

  nullify(solnVec)
  
  !!-----------------------------------------------------------------------
  !!    Update solution: small negative numbers are floored to 
  !!    gr_hypreFloor if gr_hypreUseFloor is true.
  !!    Required in the context of MGD.
  !!-----------------------------------------------------------------------  
  
  call HYPRE_SStructVectorGather(gr_hypreVecX, ierr)  
  
  var = 0
  !print*,gr_hypreSolnIsDelta;stop !this is false

  lb = 1
  call Grid_getTileIterator(itor, nodetype=blockType)
  do while(itor%isValid())

     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnVec, CENTER)

     blkLimits   = tileDesc%limits
     blkLimitsGC = tileDesc%blkLimitsGC

     level = tileDesc%level
    
     part = level - gr_hypreRefineMIN
     
     datasize(1:MDIM)=blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM)+1
     
     allocate(BoxVal(product(dataSize(1:NDIM))))
     
     if (.NOT. skipUpdate) then
        BoxVal = 0.0

        !! Use GetBoxValues more efficient then GetValues.
        call HYPRE_SStructVectorGetBoxValues(gr_hypreVecX, part,gr_hypreLower(lb,1:NDIM), &
             gr_hypreUpper(lb,1:NDIM), var, BoxVal(:), ierr)          

        do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
           do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
              do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)                      

                 ii = (k-blkLimits(LOW,KAXIS)+1)                                +  &
                      (j-blkLimits(LOW,JAXIS))*dataSize(KAXIS)                  +  &
                      (i-blkLimits(LOW,IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)

                 if (gr_hypreSolnIsDelta) then
                    solnVec(iVar,i,j,k) = BoxVal(ii) !! DEV: + ...
                 else
                    solnVec(iVar,i,j,k) = BoxVal(ii)
                 end if
                 if (gr_hypreUseFloor) then
                    if (solnVec(iVar,i,j,k) < gr_hypreFloor) solnVec(iVar,i,j,k) = gr_hypreFloor
                 end if

              end do
           end do
        end do
     end if

     ! Zero out gr_hypreVecB for the next round.
     BoxVal = 0.0
     call HYPRE_SStructVectorSetBoxValues(gr_hypreVecB, part,gr_hypreLower(lb,1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, BoxVal(:), ierr) 
     
    deallocate (BoxVal)
    call tileDesc%releaseDataPtr(solnVec, CENTER)

     lb = lb + 1
     call itor%next()

  end do
  call Grid_releaseTileIterator(itor) 
  
  call Timers_stop("gr_hypreUpdateSoln") 
  
  return
  
end subroutine gr_hypreUpdateSoln
