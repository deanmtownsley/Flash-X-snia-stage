!!****f* source/Grid/Grid_renormMassScalars
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
!!  Grid_renormMassScalars
!!
!!
!! SYNOPSIS
!!
!!  Grid_renormMassScalars(integer(IN)  :: blkLimits(2,MDIM),
!!                         real,pointer :: solnData(:,:,:,:))
!!
!! DESCRIPTION
!!
!!  Renormalize the various Mass Scalar's in groups so they sum to 1.
!!
!! ARGUMENTS
!!
!!  blkLimits - the index limits for internal zones of the block to renormalize
!!  solnData -  Pointer to the block to be renormalized
!!
!!***

subroutine Grid_renormMassScalars(blkLimits,solnData)

  implicit none

#include "constants.h"

  integer, intent(in), dimension(2,MDIM)::blkLimits
  real,pointer :: solnData(:,:,:,:)

  return
end subroutine Grid_renormMassScalars

