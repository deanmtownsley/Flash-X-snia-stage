!!****f* source/Simulation/Simulation_computeAnalytical
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
!!  Simulation_computeAnalytical
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_computeAnalytical(integer(IN) :: blockID, 
!!                                    real(IN)    :: tcurr)
!!
!!
!!
!! DESCRIPTION
!!
!!  Compute an analytical solution for a given block.
!!
!!  This is simulation-dependent, there is no general implementation.
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  tcurr   -        current time
!!
!! SIDE EFFECTS
!!
!!  The analytical solution is computed and stored in an appropriate slot
!!  (or slots) in the solution vector, UNK.
!!
!! NOTES
!! 
!!  By default this is just a stub that does not do anything.
!!
!!***

subroutine Simulation_computeAnalytical(solnData, tileDesc, tcurr)

  use Grid_tile, ONLY : Grid_tile_t
  implicit none

#include "FortranLangFeatures.fh"
  
  real,dimension(:,:,:,:),POINTER_INTENT_IN :: solnData
  type(Grid_tile_t), intent(in) :: tileDesc
  real,    intent(in) :: tcurr

end subroutine Simulation_computeAnalytical
