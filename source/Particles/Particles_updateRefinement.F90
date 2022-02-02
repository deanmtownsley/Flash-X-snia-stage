!!****f* source/Particles/Particles_updateRefinement
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
!!  Particles_updateRefinement
!!
!! SYNOPSIS
!!
!!  Particles_updateRefinement(real(inout) :: oldLocalNumBlocks)
!!
!! DESCRIPTION
!!   This routine provides a hook into the particle data structure
!!   for the Grid. It is called during update Refinement process
!!   by the Grid. The routine passes the control right back
!!   grid, with Particles specific data structures in the argument
!!   list, so that Grid can operate on them.
!!
!! ARGUMENTS
!!
!!    oldLocalNumBlocks :   number of blocks on a processor before 
!!                          refinement. 
!!
!! PARAMETERS
!!  
!!
!!***

subroutine Particles_updateRefinement(oldLocalNumBlocks)
  implicit none 
  integer,intent(INOUT) :: oldLocalNumBlocks

  return
end subroutine Particles_updateRefinement
