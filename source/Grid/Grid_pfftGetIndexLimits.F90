!!****f* source/Grid/Grid_pfftGetIndexLimits
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
!!   Grid_pfftGetIndexLimits
!!
!! SYNOPSIS
!!
!!   Grid_pfftGetIndexLimits(integer(OUT) :: configLimits(MDIM),
!!                           integer(OUT) :: phaseLimits(MDIM))
!!
!! DESCRIPTION 
!!
!!  Return the starting and ending points of the local storage
!!  in configuration and phase spaces, relative to the global
!!  domain. The configuration phase usually refers to the data
!!  at the input, and phase space to the data after forward transform
!!  
!! ARGUMENTS
!!
!!  configLimits  - endpoints in configuration space
!!  phaseLimits  - endpoints in phase space
!!
!!***
subroutine Grid_pfftGetIndexLimits(configLimits,phaseLimits)
#include "constants.h"
  implicit none 
  integer,dimension(LOW:HIGH,MDIM),intent(OUT) :: configLimits, phaseLimits
  configLimits=1
  phaseLimits=1
  return
end subroutine Grid_pfftGetIndexLimits
