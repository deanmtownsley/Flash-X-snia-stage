!!****if* source/Grid/localAPI/gr_xyzToBlock
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
!!  gr_xyzToBlock
!!
!! SYNOPSIS
!!
!!  call gr_xyzToBlock( real (IN)  :: xyz(MDIM),
!!                    integer(OUT) :: procID,
!!                    integer(OUT) :: blkID)
!!  
!! DESCRIPTION 
!!  
!!  This routine returns the identity of the block, on
!!  which the specified physical location falls. It also
!!  returns the processor ID on which the block is residing.
!!  
!!
!!
!! ARGUMENTS 
!!
!! NOTES
!!
!!  This is a stub implementation that does not do anything useful.
!!
!!  A more useful implementation is provided by PARAMESH Grid
!!  implementations. It requires PARAMESH4DEV with BITTREE.
!!***

#include "constants.h"

Subroutine gr_xyzToBlock(xyz, procID, blkID)

  implicit none

  real, dimension(MDIM),intent(IN) :: xyz
  integer, intent(OUT) :: procID
  integer, intent(OUT) :: blkID

  procID=0
  blkID=0  
End Subroutine gr_xyzToBlock
