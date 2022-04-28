!!****if* source/Particles/localAPI/pt_utAnaGetNewPosComponents
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
!!  pt_utAnaGetNewPosComponents
!!
!! SYNOPSIS
!!
!!  pt_utAnaGetNewPosComponents(real(inout)  :: particles(:,:),
!!                              integer(in) :: maxParticlesPerProc,
!!                              integer(in) :: xOutPart,
!!                              integer(in) :: yOutPart,
!!                              integer(in) :: zOutPart,
!!                              real(in)    :: t
!!                              integer(in) :: k)
!!
!! DESCRIPTION
!!
!!  Get new position components for the analytical solution
!!
!! ARGUMENTS
!!
!!  particles : the data structure containing the particles
!!              It is two dimensional real array, the first dimension
!!              represents properties associated with the data 
!!              structure, and 
!!              second dimension is index to individual elements in 
!!              the datastructure.
!! maxParticlesPerProc : the maximum count of particles allowed on 
!!                       any processor
!! xOutPart    : particle property with x component 
!! yOutPart    : particle property with y component 
!! zOutPart    : particle property with z component 
!! t        : time
!! k        : identity of a specific particle
!!
!!
!!***
subroutine pt_utAnaGetNewPosComponents(particles,maxParticlesPerProc,xOutPart,yOutPart,zOutPart,t,k)

  implicit none

#include "Simulation.h"

  integer, INTENT(in) :: maxParticlesPerProc
  real, INTENT(inout),dimension(NPART_PROPS,maxParticlesPerProc) :: particles
  integer, INTENT(in) :: xOutPart, yOutPart, zOutPart
  real, INTENT(in)           :: t
  integer, INTENT(in)           :: k


end subroutine pt_utAnaGetNewPosComponents
