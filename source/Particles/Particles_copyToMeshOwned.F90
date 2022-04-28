!!****f* source/Particles/Particles_copyToMeshOwned
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
!!  call Particles_copyToMeshOwned(particles2D)
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

subroutine Particles_copyToMeshOwned(particles2D, nprops, np, filePropNames)
 
  implicit none
  real, dimension(:,:), INTENT(in) :: particles2D
  integer, INTENT(in) :: nprops, np
  character(len=24), dimension(nprops), INTENT(in) :: filePropNames 

  return

 end subroutine Particles_copyToMeshOwned


