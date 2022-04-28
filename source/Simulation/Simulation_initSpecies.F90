!!****f* source/Simulation/Simulation_initSpecies
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
!! SYNOPSIS
!!
!!  Simulation_initSpecies()
!!
!! DESCRIPTION
!!
!!  This routine will initialize the species and species values needed for a 
!!  given setup.   The user should add the 
!!  implementation of this routine to the setups directory of a simulation 
!!  that needs to use the multispecies capabilities of the code.
!!
!!  There two general purpose implementations available in the code, one which sets standard  
!!  isotope properties for the nuclear burning source terms, and another one for the 
!!  Ionization source term.
!!
!!  This routine is called from Multispecies_init, and is called BEFORE
!!  the call to Simulation_init.  
!!
!! SEE ALSO
!!  Multispecies_init
!!  Simulation/SimulationComposition/Simulation_initSpecies
!!
!!***

subroutine Simulation_initSpecies()


implicit none
end subroutine Simulation_initSpecies
