!!****f* source/Multispecies/Multispecies_init
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
!!  Multispecies_init
!!
!! SYNOPSIS
!!
!!  Multispecies_init() 
!!                    
!!                    
!!
!! DESCRIPTION
!!  
!!  Initializes data structures of the Multispecies unit.
!!  The implementation subroutine first sets property values to UNDEFINED
!!  and then calls Simulation_initSpecies, where the user can implement
!!  setting the properties to meaningful values.
!!
!! ARGUMENTS
!!
!!  
!!  
!!  
!!
!! SEE ALSO
!!
!!   Simulation_initSpecies
!!
!!***!

subroutine Multispecies_init()

   implicit none

   

end subroutine Multispecies_init
