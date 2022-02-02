!!****f* source/Particles/Particles_accumCount
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
!!    Particles_accumCount
!!
!! SYNOPSIS
!!
!!    Particles_accumCount(integer, intent(IN) :: var)
!!
!! DESCRIPTION
!!    
!!   This routine is to be used in the refinement based upon the number of
!!   particles in a block. It adds a weight to the cell of the grid in the
!!   specified grid variable if a particle is found to be contained in the 
!!   cell 
!!
!! ARGUMENTS
!!
!!  var:        the grid variable to add the weights if particle found in the cell
!!
!!
!!***


subroutine Particles_accumCount(var)

  implicit none

  integer, intent(IN) :: var
!----------------------------------------------------------------------

  
  return
  
  !----------------------------------------------------------------------
  
end subroutine Particles_accumCount


