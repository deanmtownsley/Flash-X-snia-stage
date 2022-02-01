!!****f* source/Simulation/Simulation_getVarnameType
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
!!  Simulation_getVarnameType
!!
!! SYNOPSIS
!!
!!  Simulation_getVarnameType(integer, intent(in)  :: varname,
!!                            integer, intent(out)  :: vartype)
!!
!! DESCRIPTION
!!  Given a variable name (e.g., DENS_VAR), this subroutine
!!  returns the type in vartype as one of GENERIC, PER_VOLUME,
!!  or PER_MASS.
!!
!!  If the simulation is not configured correctly, VARTYPE_ERROR
!!  may be returned instead.
!!
!! ARGUMENTS
!!
!!   varname : name of variable
!!
!!   vartype : type of variable varname.
!!
!! NOTES
!!  An implementation of this subroutine is usually automatically
!!  generated at setup time and should not be edited.
!!
!!***


subroutine Simulation_getVarnameType(varname,vartype)
implicit none 

#include "constants.h"
#include "Simulation.h"

   integer, intent(out) :: vartype
   integer, intent(in) :: varname

   vartype = VARTYPE_ERROR  ! By default vartype is an ERROR

end subroutine Simulation_getVarnameType


! Stub file for getVarnameType
! given a variable name e.g. DENS_VAR
! returns the type it belongs to VARTYPE_{GENERIC,PER_VOLUME,PER_MASS}
