!!****if* source/Grid/GridSolvers/HYPRE/gr_hypreFinalize
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
!!  gr_hypreFinalize
!!
!!
!! SYNOPSIS
!!
!!  call gr_hypreFinalize ()
!!
!! Description
!!  
!!  Release HYPRE objects and other temporary storage from memory.
!!
!! ARGUMENTS
!!
!!  none  
!!
!! PARAMETERS
!!
!!***

subroutine gr_hypreFinalize()

  use gr_hypreLocalInterface, only: gr_hypreDestroyGrid, gr_hypreDestroySolver  
  implicit none
  
  !! Destroy the HYPRE grid object.
  call gr_hypreDestroyGrid()
  call gr_hypreDestroySolver()
  
end subroutine gr_hypreFinalize
