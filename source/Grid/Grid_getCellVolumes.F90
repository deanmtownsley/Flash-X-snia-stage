!!****f* source/Grid/Grid_getCellVolumes
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

