!!****f* source/Grid/Grid_conserveToPrimitive
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
!!  Grid_conserveToPrimitive
!!
!!
!! SYNOPSIS
!!
!!  call Grid_conserveToPrimitive(integer(in) :: blkList(count),
!!                           integer(in) :: count,
!!                           logical(in) :: allCells,
!!                           logical(in) :: force)
!!
!!
!! DESCRIPTION
!!
!!  Calls gr_conserveToPrimitive
!!
!!
!! ARGUMENTS
!! 
!!   blkList - integer list of blocks to be operated on
!!
!!   count - number of blocks in the blkList
!!
!!   allCells - act on all cells, including guardcells, if .TRUE.,
!!              otherwise only modify interior cells.
!!
!!   force - whether to force conversion
!!
!! NOTES
!!
!!  DEV: Currently only implemented for Paramesh4!
!!***

subroutine Grid_conserveToPrimitive(blkList,count,allCells,force)
  implicit none
  integer,intent(IN) :: count
  integer,dimension(count),intent(IN) :: blkList 
  logical,intent(IN) :: allCells, force
end subroutine Grid_conserveToPrimitive
