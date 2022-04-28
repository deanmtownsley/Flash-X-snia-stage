!!****f* source/Grid/Grid_getCellFaceAreas
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
!!  Grid_getCellFaceAreas
!!
!! SYNOPSIS
!!  call Grid_getCellFaceAreas(integer(IN) :: axis,
!!                             integer(IN) :: level,
!!                             integer(IN) :: lo(1:MDIM),
!!                             integer(IN) :: hi(1:MDIM),
!!       real(OUT),dimension(lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)) :: areas)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!   axis - specifies the face-centered index space 
!!   level - refinement level.
!!           This is 1-based, i.e., the root level is numbered 1.
!!   lo,hi - the lower-left and upper-right points in the global face-centered
!!           index space that specify the region of faces whose areas are 
!!           requested
!!   areas - the array in which the requested area values will be stored
!!
!!***

#include "constants.h"

subroutine Grid_getCellFaceAreas(axis, level, lo, hi, areas)
   implicit none
   integer, intent(IN)  :: axis
   integer, intent(IN)  :: level
   integer, intent(IN)  :: lo(1:MDIM)
   integer, intent(IN)  :: hi(1:MDIM)
   real,    intent(OUT) :: areas(lo(IAXIS):hi(IAXIS), &
                                 lo(JAXIS):hi(JAXIS), &
                                 lo(KAXIS):hi(KAXIS))

   areas(:,:,:) = 0.0
end subroutine Grid_getCellFaceAreas

