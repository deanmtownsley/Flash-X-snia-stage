!!****f* source/Particles/Particles_copyFromMeshOwned
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
!!  Grid_copyFromMeshOwned
!!
!! SYNOPSIS
!!
!!  call Particles_copyFromMeshOwned(np, particleType, meshOwnedPart, particles2D)
!!                    
!!  
!! DESCRIPTION 
!!  
!! Get position, velocity, id, cpu from amrex particle and 
!! populate particles2D for writing into hdf5 file later
!!
!! ARGUMENTS 
!!
!!  np : number of particles
!!
!!  meshOwnedPart : particles in mesh-owned format
!!
!!  particles2D : 2D particles with property values
!!
!!  particleType : Type of particle : passive, active, star etc.
!!
!!***
!===============================================================================

subroutine Particles_copyFromMeshOwned(np, particleType, meshOwnedPart, particles2D)

  use Particles_dataTypes, ONLY : Particles_meshParticle_t
  
  implicit none

  integer, INTENT(in) :: np, particleType
  type(Particles_meshParticle_t), INTENT(in) :: meshOwnedPart(:)
  real, dimension(:,:), INTENT(out)          :: particles2D

  particles2D( :, :) = 0.0
  return

 end subroutine Particles_copyFromMeshOwned


