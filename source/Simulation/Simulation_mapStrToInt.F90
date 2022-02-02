!!****f* source/Simulation/Simulation_mapStrToInt
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
!!  Simulation_mapStrToInt
!!
!!
!! SYNOPSIS
!!  Simulation_mapStrToInt(character, len=*, intent(IN) :: str,
!!                         integer, intent(OUT)         :: key, 
!!                         integer, intent(IN)          :: map)
!!
!!
!! DESCRIPTION
!!
!!  This routine is created by the setup script, and should never be edited.
!!  This routine maps a string described in the Config file to an integer
!!  index described in the Simulation.h file.  The integer can represent a variable,
!!  a species, a flux variable, or a particle property.  
!!
!!  For example, the Config file might say:
!!  VARIABLE velx
!!  At setup time, the Simulation.h file would be created to look like
!!  #define VELX_VAR 6
!!  The result of "call Simulation_mapStrToInt(keyout,'velx',MAPBLOCK_UNK)" would
!!  be keyout=6.    
!!
!! ARGUMENTS
!! 
!!  str   --  string input
!!  key   --  returned integer index
!!  map   --  variable indicating the type of data structure within Simulation.h.  Valid values are
!!            MAPBLOCK_UNK   for variables and species (NAME_VAR or NAME_SPEC in Simulation.h)
!!            MAPBLOCK_FLUX  for flux variables (NAME_FLUX in Simulation.h)          
!!            MAPBLOCK_PART  for particle properties (NAME_PART_PROP in Simulation.h)
!!
!!***

subroutine Simulation_mapStrToInt(str,key,map)
implicit none 

#include "constants.h"

   character(len=*), intent(in) :: str
   integer, intent(out) :: key 
   integer, intent(IN) :: map

   key = NONEXISTENT

end subroutine Simulation_mapStrToInt

