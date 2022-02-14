!!****ih* source/Grid/GridSolvers/HYPRE/gr_hypreLocalInterface
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
!! This is a module file for the HYPRE solver in FLASH that defines
!! additional interfaces private to the GridSolvers/HYPRE implementations.
!!
!! NOTES
!!
!!  Adding explicit interfaces here is being tried as an alternative to
!!  writiting executable FORTRAN wrappers for each additional HYPRE routine
!!  we want to call from FLASH.
!!
!! SEE ALSO
!!
!!  gr_hypreF90CAdapters.c
!!***

#include "constants.h"

Module gr_hypreLocalInterface

#if 0
  ! Maybe one day... for now, not all compilers support this.
  interface
     subroutine hypre_sstructinnerprod(fx, fy, fresult, ierr) &
          bind(C,name="HYPRE_SStructInnerProd")
       implicit none
       integer*8,intent(IN)  :: fx, fy
       real,     intent(OUT) :: fresult
       integer,  intent(OUT) :: ierr
    end subroutine hypre_sstructinnerprod

     integer function hypre_pcggetconverged(solver, converged) &
          bind(C,name="HYPRE_PCGGetConverged")
       implicit none
       integer*8,VALUE :: solver ! intent(IN)
       integer,intent(OUT) :: converged
     end function hypre_pcggetconverged

     subroutine hypre_describeerror(hypreErr, description) &
          bind(C,name="HYPRE_DescribeError")
       implicit none
       integer,VALUE :: hypreErr ! intent(IN)
       character(len=1),dimension(MAX_STRING_LENGTH),intent(OUT) :: description
     end subroutine hypre_describeerror

  end interface

#else
  ! fallback... wrapper implementations are in gr_hypreF90CAdapters.c

  interface
     subroutine hypre_sstructinnerprod(fx, fy, fresult, ierr)
       implicit none
       integer*8,intent(IN)  :: fx, fy
       real,     intent(OUT) :: fresult
       integer,  intent(OUT) :: ierr
    end subroutine hypre_sstructinnerprod

     subroutine hypre_pcggetconverged(solver, converged, ierr)
       implicit none
       integer*8,intent(IN) :: solver
       integer,intent(OUT) :: converged
       integer,  intent(OUT) :: ierr
     end subroutine hypre_pcggetconverged

     subroutine hypre_describeerror(hypreErr, description)
       implicit none
       integer,intent(IN) :: hypreErr
       character(len=MAX_STRING_LENGTH),intent(OUT) :: description
     end subroutine hypre_describeerror

  end interface

#endif

#include "FortranLangFeatures.fh"
#include "Simulation.h"  

  interface
     subroutine gr_hypreApplyBcToFace(blkLimits,blkLimitsGC,part,var,bcType,direction, &
                                 bcValue, del, Lower, tileDesc)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       integer, intent(IN) :: blkLimits (2,MDIM) 
       integer, intent(IN) :: blkLimitsGC (2,MDIM)
       integer, intent(IN) :: part
       integer, intent(IN) :: var
       integer, intent(IN) :: bcType
       integer, intent(IN) :: direction
       real,    intent(IN) :: bcValue(2)
       real,    intent(IN) :: del
       integer, intent(IN) :: Lower(MDIM)
       type(Grid_tile_t), intent(IN) :: tileDesc       
     end subroutine gr_hypreApplyBcToFace
     subroutine gr_hypreApplyBcToFace_VecB(blkLimits,blkLimitsGC,part,var,bcType,direction, &
                                 bcValue, del, Lower, tileDesc)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       integer, intent(IN) :: blkLimits (2,MDIM) 
       integer, intent(IN) :: blkLimitsGC (2,MDIM)
       integer, intent(IN) :: part
       integer, intent(IN) :: var
       integer, intent(IN) :: bcType
       integer, intent(IN) :: direction
       real,    intent(IN) :: bcValue(2)
       real,    intent(IN) :: del
       integer, intent(IN) :: Lower(MDIM)
       type(Grid_tile_t), intent(IN) :: tileDesc       
     end subroutine gr_hypreApplyBcToFace_VecB
  end interface

  interface
     subroutine gr_hypreGridStatus (blockCount, blockType, nvars)
       implicit none  
       integer, intent(IN):: blockCount
       integer, intent(IN):: blockType
       integer, intent(IN),OPTIONAL :: nvars
     end subroutine gr_hypreGridStatus
  end interface
  
  interface
     subroutine gr_hypreSetupGrid (blockCount, blockType, nvars)
       implicit none 
       integer,                      intent(IN) :: blockCount
       integer,                      intent(IN) :: blockType
       integer,OPTIONAL,             intent(IN) :: nvars
     end subroutine gr_hypreSetupGrid
  end interface
 
  interface
     subroutine gr_hypreComputeB (iVar, blockCount, blockType, bcTypes, bcValues)
       integer, intent(IN) :: iVar
       integer, intent(IN) :: blockCount
       integer, intent(IN) :: blockType
       integer, intent(IN) :: bcTypes(6)
       real,    intent(IN) :: bcValues(2,6)  
     end subroutine gr_hypreComputeB
  end interface
  
  interface 
     subroutine gr_hypreCreateMatrix(blockCount, blockType, bcTypes, bcValues)
       integer, intent(IN) :: bcTypes(6)
       real,    intent(IN) :: bcValues(2,6)
       integer, intent(IN) :: blockCount
       integer, intent(IN) :: blockType
     end subroutine gr_hypreCreateMatrix
  end interface

  interface
     subroutine gr_hypreSetIniGuess (iVar, blockCount, blockType)
       integer, intent(IN) :: iVar
       integer, intent(IN) :: blockCount
       integer, intent(IN) :: blockType
     end subroutine gr_hypreSetIniGuess
  end interface

  interface
     subroutine gr_hypreUpdateSoln (iVar, blockCount, blockType)
       integer, intent(IN) :: iVar
       integer, intent(IN) :: blockCount
       integer, intent(IN) :: blockType
     end subroutine gr_hypreUpdateSoln
  end interface

  interface
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
     end subroutine gr_hypreAddGraph
  end interface

  interface 
     subroutine gr_hypreDestroyGrid()
     end subroutine
  end interface
 
  interface 
     subroutine gr_hypreSetupSolver()
     end subroutine
  end interface
 
  interface 
     subroutine gr_hypreDestroySolver()
     end subroutine
  end interface
 
  interface 
     subroutine gr_hypreSolve()
     end subroutine
  end interface

  interface
     subroutine gr_hypreLevel()
     end subroutine
  end interface

end Module gr_hypreLocalInterface
