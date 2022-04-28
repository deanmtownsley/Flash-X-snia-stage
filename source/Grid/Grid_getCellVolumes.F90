!!****f* source/Grid/Grid_getCellVolumes
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
!!  Grid_getCellVolumes
!!
!! SYNOPSIS
!!  call Grid_getCellVolumes(integer(IN) :: level,
!!                           integer(IN) :: lo(1:MDIM),
!!                           integer(IN) :: hi(1:MDIM),
!!       real(OUT),dimension(lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS))::volumes)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!   level - refinement level.
!!           This is 1-based, i.e., the root level is numbered 1.
!!   lo    - the lower-left point in the global cell-centered index
!!           space that specifies, togehter with hi,  the region of
!!           cells whose volumes are requested.
!!   hi    - the upper-right point in the global cell-centered index
!!           space that specifies, together with lo,  the region of
!!           cells whose volumes are requested.
!!   volumes - the array in which the requested volume values will be stored
!!
!!***

#include "constants.h"

subroutine Grid_getCellVolumes(level, lo, hi, volumes)
   implicit none
   integer, intent(IN)  :: level
   integer, intent(IN)  :: lo(1:MDIM)
   integer, intent(IN)  :: hi(1:MDIM)
   real,    intent(OUT) :: volumes(lo(IAXIS):hi(IAXIS), &
                                   lo(JAXIS):hi(JAXIS), &
                                   lo(KAXIS):hi(KAXIS))

   volumes(:,:,:) = 0.0
end subroutine Grid_getCellVolumes

