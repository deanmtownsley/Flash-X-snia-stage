!!****if* source/Grid/GridSolvers/HYPRE/paramesh/gr_hypreLevel
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
!!    gr_hypreLevel
!!
!!  SYNOPSIS
!!   
!!      call gr_hypreLevel()
!!
!!  DESCRIPTION 
!!  
!! NOTES:
!!
!!   Uses HYPRE library.  
!!
!!***
subroutine gr_hypreLevel ()

   use tree, only: lrefine
   use gr_hypreData, only: gr_hypreRefineMAX, gr_hypreRefineMIN

   implicit none

   gr_hypreRefineMAX = max(maxval(lrefine), gr_hypreRefineMAX)
   gr_hypreRefineMIN = 1

end subroutine gr_hypreLevel
