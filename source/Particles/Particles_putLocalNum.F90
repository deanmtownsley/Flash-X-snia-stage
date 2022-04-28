!!****f* source/Particles/Particles_putLocalNum
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
!!  Particles_putLocalNum
!!
!! SYNOPSIS
!!
!!  Particles_putLocalNum(integer(IN)  :: localNumParticles)
!!               
!!  
!! DESCRIPTION 
!!
!!  Sets the local number of particles on a given proc
!!  
!!
!! ARGUMENTS
!!
!!  localNumParticles:  the local number of particles
!!
!! NOTES
!!
!!
!!***

subroutine Particles_putLocalNum(localNumParticles)

implicit none
  integer, intent(in)  :: localNumParticles

  return
end subroutine Particles_putLocalNum
