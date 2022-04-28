!!****f* source/Particles/Particles_updateAttributes
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
!!  Particles_updateAttributes
!!
!! SYNOPSIS
!!
!!  Particles_updateAttributes()
!!
!!  Particle attribute advancement routines.  Use this routine to 
!!  update user-define particle attributes beyond the usual
!!  particle attributes of position, velocity, block, tag, and mass.
!!  It is usually called from  IO_output. It makes sure that guardcells
!!  are current and, if necessary, that Eos has been applied to them
!!  before particles attributes are calculated from the grid.
!!   
!!  The attributes can be specified at runtime 
!!  using runtime parameters particle_attribute_1, particle_attribute_2
!!  etc.
!!
!! ARGUMENTS
!!  
!! NOTES
!!
!! The map between particle property and a mesh variable on which the
!! property is dependent has to be specified in the Config file of 
!! the Simulation directory for this routine to work right. Please
!! see the Config file of IsentropicVortex setup for an example, and
!! also see the Setup chapter of the User's Guide.
!! 
!!
!!  
!!
!!
!!***


subroutine Particles_updateAttributes()


  implicit none


end subroutine Particles_updateAttributes

