!!****f* source/Grid/Grid_getNeighProcList
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
!!  Grid_getNeighProcList
!!
!! SYNOPSIS
!!
!!  call Grid_getNeighProcList(logical(IN)             :: includeMyProc, 
!!                             integer(INOUT), pointer :: neighProcList(:),
!!                             integer(OUT)            :: numNeigh)
!!
!! DESCRIPTION 
!!
!!  Creates a pointer array containing the neighboring processor IDs of all
!!  LEAF blocks on this processor.
!!
!! ARGUMENTS 
!!  
!!
!!  includeMyProc - Whether the array should include my processor ID.
!!
!!  neighProcList - The processor IDs of all neighboring LEAF blocks. 
!!
!!  numNeigh      - The number of entries in neighProcList.
!!
!!
!! NOTES
!!
!!  Currently only implemented for Paramesh4 Grid implementations.
!!
!!  It is the users resposibility to
!!    1. deallocate neighProcList
!!    2. obtain a new neighProcList when the mesh changes
!!
!!***

subroutine Grid_getNeighProcList(includeMyProc, neighProcList, numNeigh)
  implicit none
  logical, intent(IN) :: includeMyProc
  integer, dimension(:), pointer :: neighProcList
  integer, intent(OUT) :: numNeigh
  numNeigh = 0
end subroutine Grid_getNeighProcList
