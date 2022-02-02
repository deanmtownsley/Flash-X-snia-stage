!!****if* source/Grid/GridSolvers/HYPRE/gr_hypreAddGraph
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
!!    gr_hypreAddGraph
!!
!!  SYNOPSIS
!!   
!!      call gr_hypreAddGraph()
!!
!!  DESCRIPTION 
!!    This routine handles addition of Graphs (fine-coarse boundaries) across a
!!    a particular face.
!! 
!!
!!  
!! NOTES:
!!
!!   Uses HYPRE library.  
!!
!!***
      
subroutine gr_hypreAddGraph (blockHandle, blockID, blkPartNo, direction, datasize, &
                             CornerID, blkStride, firstHypreVar, numVars)
  
      integer, intent(IN) :: blockHandle
      integer, intent(IN) :: direction
      integer, intent(IN) :: blockID
      integer, intent(IN) :: blkPartNo
      integer, intent(IN) :: datasize(MDIM)
      integer, intent(IN) :: CornerID(MDIM)
      integer, intent(IN) :: blkStride(MDIM)
      integer, intent(IN),OPTIONAL :: firstHypreVar, numVars

      return
end subroutine gr_hypreAddGraph
