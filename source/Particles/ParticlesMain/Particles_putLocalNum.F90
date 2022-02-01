!!****if* source/Particles/ParticlesMain/Particles_putLocalNum
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
!!  
!! DESCRIPTION 
!!
!!  Sets the number of particles local to the current MPI task
!!  
!!
!! ARGUMENTS
!!
!!  localNumParticles - input value of local number of particles
!!
!! NOTES
!!
!!
!!***

subroutine Particles_putLocalNum(localNumParticles)

  use Particles_data, ONLY : pt_numLocal

implicit none
  integer, intent(in)  :: localNumParticles

  pt_numLocal = localNumParticles


  return
end subroutine Particles_putLocalNum
