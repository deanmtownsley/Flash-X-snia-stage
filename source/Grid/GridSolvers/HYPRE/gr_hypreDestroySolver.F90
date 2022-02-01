!!****if* source/Grid/GridSolvers/HYPRE/gr_hypreDestroySolver
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
!! gr_hypreDestroySolver
!!
!!  SYNOPSIS
!!
!!  call gr_hypreDestroySolver ()
!!
!!
!!  DESCRIPTION 
!! This routine destroys the HYPRE solver object 
!! with it's associated PC. This routine is called only once.
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
!!
!!***

!!REORDER(4): solnVec

subroutine gr_hypreDestroySolver ()
  
  use gr_hypreData, ONLY : gr_hypreSolver, gr_hyprePC, &
                           gr_hyprePCType, gr_hypreSolverType     

  use Timers_interface, ONLY : Timers_start, Timers_stop  

  implicit none
  
#include "Simulation.h"
#include "constants.h"   
#include "HYPREf.h" 

  integer :: ierr

  call Timers_start("gr_hypreDestroySolver")    

  select case (gr_hypreSolverType)
  case (HYPRE_PCG)
     call HYPRE_ParCSRPCGDestroy(gr_hypreSolver, ierr)        
  case (HYPRE_BICGSTAB)
     call HYPRE_parCSRBiCGSTABDestroy(gr_hypreSolver, ierr)    
  case (HYPRE_AMG)
     call HYPRE_BoomerAMGDestroy(gr_hypreSolver, ierr)
  case (HYPRE_GMRES)
     call HYPRE_parCSRGMRESDestroy(gr_hypreSolver, ierr)    
  case (HYPRE_SPLIT)
     call HYPRE_SStructSplitDestroy(gr_hypreSolver, ierr)
  case (HYPRE_HYBRID)     
     call HYPRE_ParCSRHybridDestroy(gr_hypreSolver, ierr)
  end select
  
  select case (gr_hyprePCType)
  case (HYPRE_AMG)
     call HYPRE_BoomerAMGDestroy(gr_hyprePC, ierr)     
  case (HYPRE_ILU)
     call HYPRE_EuclidDestroy(gr_hyprePC, ierr)
  case (HYPRE_PARASAILS)
     call HYPRE_ParaSailsDestroy(gr_hyprePC, ierr)     
  end select
  
  call Timers_stop("gr_hypreDestroySolver")    
  
  
end subroutine gr_hypreDestroySolver
