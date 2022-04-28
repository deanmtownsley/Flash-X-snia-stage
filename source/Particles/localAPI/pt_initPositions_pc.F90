!!****if* source/Particles/localAPI/pt_initPositions_pc
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
!!    pt_initPositions
!!
!! SYNOPSIS
!!
!!    call pt_initPositions(integer(in)  :: blockID,
!!                          logical(out) :: success)
!!
!! DESCRIPTION
!!
!!  Initializes particle locations for one block in the grid.
!!
!!  The initialization of locations may fail, if adding the
!!  number of particles requested (by the choice of initialization
!!  method and runtime parameters) would exceed the maximum allowed
!!  number of particles allowed in any MPI task given by runtime
!!  parameter pr_maxPerProc. In that case, an impementation of this
!!  interface may:
!!    o abort the simulation (by calling Driver_abortFlash); or
!!    o return with success set to FALSE.
!!  If returning with success=FALSE, pt_numLocal may have been
!!  reset to 0.
!!
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!!***


subroutine pt_initPositions_pc (iSpecies,tileDesc,success)

  use Grid_tile,        ONLY : Grid_tile_t
  implicit none

  integer, intent(in) :: iSpecies
  type(Grid_tile_t), INTENT(in) :: tileDesc
  logical,intent(OUT) :: success

  success = .FALSE.
  return

!----------------------------------------------------------------------
  
end subroutine pt_initPositions_pc


