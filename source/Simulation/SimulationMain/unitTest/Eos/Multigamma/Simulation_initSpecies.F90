!!****if* source/Simulation/SimulationMain/unitTest/Eos/Multigamma/Simulation_initSpecies
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
!!  Simulation_initSpecies
!!
!!
!! SYNOPSIS
!!  Simulation_initSpecies()
!!
!! DESCRIPTION
!!
!!  This routine will initialize the simulation with the species 
!!  in the Config files. This particular implementation has air
!!  and SF6, and they have two different gamma values for the 
!!  ideal gas gamma law. 
!!
!!
!!
!!***

subroutine Simulation_initSpecies()
  use Multispecies_interface, ONLY : Multispecies_setProperty

implicit none
#include "Simulation.h"
#include "Multispecies.h"


  call Multispecies_setProperty(HE4_SPEC, A, 4.)
  call Multispecies_setProperty(HE4_SPEC, Z, 2.)
  call Multispecies_setProperty(HE4_SPEC, GAMMA, 1.66)

  call Multispecies_setProperty(NI56_SPEC, A, 56.)
  call Multispecies_setProperty(NI56_SPEC, Z, 28.)
  call Multispecies_setProperty(NI56_SPEC, GAMMA, 1.4)

end subroutine Simulation_initSpecies
