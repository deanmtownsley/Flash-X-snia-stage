!!****if* source/Grid/GridMain/AMR/Amrex/Grid_primitiveToConserve
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
!!  Grid_primitiveToConserve
!!
!!
!! SYNOPSIS
!!  Grid_primitiveToConserve(integer(in) :: blkList(count),
!!                           integer(in) :: count,
!!                           logical(in) :: force)
!!
!!
!! DESCRIPTION
!!  Calls gr_primitiveToConserve
!!
!! ARGUMENTS
!!   blkList - integer list of blocks to be operated on
!!   count - number of blocks in the blkList
!!   force - whether to force conversion
!!
!!***

subroutine Grid_primitiveToConserve(blkList, count, force)
  use Driver_interface, ONLY : Driver_abort

  implicit none

  integer, intent(IN) :: count
  integer, intent(IN) :: blkList(count)
  logical, intent(IN) :: force

  ! DEV: TODO This needs to be rethought or modernized to work with iterators
  call Driver_abort("[Grid_primitiveToConserve] Not implemented for AMReX")
end subroutine Grid_primitiveToConserve

