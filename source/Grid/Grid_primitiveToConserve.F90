!!****f* source/Grid/Grid_primitiveToConserve
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
!!  Grid_primitiveToConserve
!!
!!
!! SYNOPSIS
!!
!!  call Grid_primitiveToConserve(integer(in) :: blkList(count),
!!                           integer(in) :: count,
!!                           logical(in) :: force)
!!
!!
!! DESCRIPTION
!!
!!  Calls gr_primitiveToConserve
!!
!!
!! ARGUMENTS
!! 
!!   blkList - integer list of blocks to be operated on
!!
!!   count - number of blocks in the blkList
!!
!!   force - whether to force conversion
!!
!! NOTES
!!
!!  DEV: Currently only implemented for Paramesh4!
!!***

subroutine Grid_primitiveToConserve(blkList,count,force)
  implicit none
  integer,intent(IN) :: count
  integer,dimension(count),intent(IN) :: blkList 
  logical,intent(IN) :: force
end subroutine Grid_primitiveToConserve
