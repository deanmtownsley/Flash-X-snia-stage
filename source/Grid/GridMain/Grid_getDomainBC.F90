!!****if* source/Grid/GridMain/Grid_getDomainBC
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
!!  Grid_getDomainBC
!!
!! SYNOPSIS
!!
!!  call Grid_getDomainBC(integer(OUT) :: boundary(LOW:HIGH,MDIM))
!!                    
!! DESCRIPTION 
!!  Returns the boundary condition for each face of the domain
!!  
!!  The boundary conditions are defined in the header file
!!  constants.h, ie. OUTFLOW, REFLECTING, PERIODIC etc
!!   
!! ARGUMENTS 
!!
!!  boundary   - array returned holding boundary conditions, except when periodic,
!!            if any of the faces of the block are on a physical boundary. 
!!                
!!            the first index of the array can take on values LOW or
!!            HIGH, and the second index can be IAXIS, JAXIS or KAXIS
!!
!! NOTES
!!
!!   The #define constants LOW, HIGH, IAXIS, JAXIS and KAXIS
!!   are defined in constants.h and are
!!   meant to ease the readability of the code.       
!!   instead of boundary(2,3) = PERIODIC the code reads
!!   boundary(HIGH, KAXIS) = PERIODIC
!!
!!***


subroutine Grid_getDomainBC(boundary)

  use Grid_data, ONLY : gr_domainBC
#include "constants.h"
  implicit none
  integer, dimension(LOW:HIGH,MDIM),intent(out):: boundary

  boundary(LOW:HIGH,1:MDIM) = gr_domainBC(LOW:HIGH,1:MDIM)
 
  return
end subroutine Grid_getDomainBC





