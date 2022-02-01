!!****f* source/Simulation/Simulation_mapIntToStr
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
!!  Simulation_mapIntToStr
!!
!!
!! SYNOPSIS
!!  Simulation_mapIntToStr(integer, intent(IN)             :: key, 
!!                         character, len=*, intent(INOUT) :: str,
!!                         integer, intent(IN)             :: map)
!!
!!
!! DESCRIPTION
!!
!!  This routine is created by the setup script, and should never be edited.
!!  This routine maps an index described in the Simulation.h file to a string
!!  described in the Config file.  The integer can represent a variable,
!!  a species, a flux variable, or a particle property.  
!!
!!  For example, the Config file might say:
!!  VARIABLE velx
!!  At setup time, the Simulation.h file would be created to look like
!!  #define VELX_VAR 6
!!  The result of "call Simulation_mapIntToStr(VELX_VAR,result,MAPBLOCK_UNK)" would
!!  be result="velx"    
!!
!! ARGUMENTS
!! 
!!  key   --  integer index
!!  str   --  returned string
!!  map   --  variable indicating the type of data structure within Simulation.h.  Valid values are
!!            MAPBLOCK_UNK   for variables and species (NAME_VAR or NAME_SPEC in Simulation.h)
!!            MAPBLOCK_FLUX  for flux variables (NAME_FLUX in Simulation.h)          
!!            MAPBLOCK_PART  for particle properties (NAME_PART_PROP in Simulation.h)
!!
!!***

subroutine Simulation_mapIntToStr(key, str, block)
implicit none 

#include "constants.h"

   integer, intent(in) :: key, block
   character(len=*), intent(inout) :: str
  
   str = "ERROR"

end subroutine Simulation_mapIntToStr


