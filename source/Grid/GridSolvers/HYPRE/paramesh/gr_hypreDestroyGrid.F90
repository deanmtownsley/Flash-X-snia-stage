!!****if* source/Grid/GridSolvers/HYPRE/paramesh/gr_hypreDestroyGrid
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
!! gr_hypreDestroyGrid
!!
!!  SYNOPSIS
!!
!!  call gr_hypreDestroyGrid()
!!
!!
!!  DESCRIPTION 
!! This routine destroys HYPRE Grid. 
!! 
!! 
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  
!! Uses HYPRE library.
!!
!!***

!!REORDER(4): solnVec


subroutine gr_hypreDestroyGrid()

  use gr_hypreData,     ONLY : gr_hypreLower, gr_hypreUpper, &
                           gr_hypreGrid, gr_hypreStencils, &
                               gr_hypreGraph, gr_hypreMatA, gr_hypreVecB, gr_hypreVecX, &
                               gr_hypreGridIsSetUp, &
                               gr_hypreNeghLevels, gr_hypreSurrBlkSum
  
  implicit none
  
#include "Simulation.h"  
#include "HYPREf.h"     
  
  integer :: ierr, s
  
  if (gr_hypreGridIsSetUp) then
     deallocate (gr_hypreLower)
     deallocate (gr_hypreUpper)
     deallocate (gr_hypreNeghLevels)     
     deallocate (gr_hypreSurrBlkSum)
     
     !! HYPRE clean up
     call HYPRE_SStructGridDestroy(gr_hypreGrid, ierr)
     if (allocated(gr_hypreStencils)) then
        do s=0,size(gr_hypreStencils)-1
           call HYPRE_SStructStencilDestroy(gr_hypreStencils(s), ierr)
        end do
        deallocate(gr_hypreStencils)
     end if
     call HYPRE_SStructGraphDestroy(gr_hypreGraph, ierr)  
     call HYPRE_SStructMatrixDestroy(gr_hypreMatA, ierr)  
     call HYPRE_SStructVectorDestroy(gr_hypreVecB, ierr)
     call HYPRE_SStructVectorDestroy(gr_hypreVecX, ierr)     
     
     gr_hypreGridIsSetUp = .FALSE.
  end if
  
  
end subroutine gr_hypreDestroyGrid
