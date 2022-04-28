!!****if* source/Grid/localAPI/gr_bcApplyToAllBlks
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
!!  gr_bcApplyToAllBlks
!!
!! SYNOPSIS
!!
!!  gr_bcApplyToAllBlks(integer(IN) :: axis,
!!                      logical(IN) :: isWork)
!!  
!! DESCRIPTION 
!!
!!  This routine is a wrapper around gr_bcApplyToOneFace, and is used by UG and PM2.
!!  It calls gr_bcApplyToOneFace one each for lowerface and upperface, and repeats
!!  the process for all blocks in the grid.
!!  
!! 
!! ARGUMENTS
!!  
!!    axis           - the direction for applying BC, one of IAXIS, JAXIS, or KAXIS
!!    isWork         - is always false for UG. In PM2, if true, indicates that
!!                     the boundary conditions should be applied to work array
!!
!! NOTES
!!  A specific direction is required in axis - no ALLDIR at this point.
!!
!!***

subroutine gr_bcApplyToAllBlks(axis,isWork)
  implicit none
  
  integer, intent(in) :: axis
  logical, intent(in) :: isWork
  
end subroutine gr_bcApplyToAllBlks
