!!****if* source/Particles/ParticlesMain/passive/unitTest/pt_utData
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
!!  pt_utData
!!
!! SYNOPSIS
!!  use pt_utData
!!
!! DESCRIPTION
!!
!!  Store the simulation data for unitTesting of Particles
!!   
!! ARGUMENTS
!!
!! PARAMETERS   
!!      Described in the Config file
!!
!!
!!***

module pt_utData

  implicit none
  
  real, save    :: pt_utA0, pt_utA1, pt_utInitialSimTime

  logical,save :: pt_utAnalyticParticlePositions, pt_utFakeMapMeshToParticles

end module pt_utData
