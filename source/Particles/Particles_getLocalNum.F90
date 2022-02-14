!!****f* source/Particles/Particles_getLocalNum
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
!!  Particles_getLocalNum
!!
!! SYNOPSIS
!!
!!  Particles_getLocalNum(integer(IN)  :: blockID,
!!                        integer(OUT) :: localNumParticles)
!!               
!!  
!! DESCRIPTION 
!!
!!  Returns the local number of particles on a given block
!!     Needed for input/output routines
!!  
!! ARGUMENTS
!!
!!  blockID - block number !!DEV: IGNORED!
!!  localNumParticles - returned value of local particles
!!
!! NOTES
!!
!!
!!***

subroutine Particles_getLocalNum(blockID, localNumParticles)

implicit none
  integer, intent(in) :: blockID
  integer, intent(out)  :: localNumParticles

  !return 0 particles for stub level
  localNumParticles = 0

  return
end subroutine Particles_getLocalNum
