!!****f* source/Grid/Grid_receiveInputData
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
!!  Grid_receiveInputData
!!
!! SYNOPSIS
!!
!!  call Grid_receiveInputData(integer(IN) :: localNumBlocks,
!!                             integer(IN) :: alnblocks,
!!                             integer(IN) :: xx)
!!
!! DESCRIPTION 
!!
!!  Initializes grid arrays from arrays read by the I/O unit.
!!
!! ARGUMENTS  
!!
!!  localNumBlocks : the number of blocks on my processor.
!!
!!  alnblocks : the approximate number of local blocks on each
!!              processor if we give each processor an equal
!!              number of blocks.  Calculated from
!!              int(globalNumBlocks/meshNumProcs) + 1.
!!
!!  xx : an integer representing a cutoff point.  Processors
!!       less than this value are assigned alnblocks blocks and
!!       processors greater than or equal to this value are
!!       assigned lnblocks-1 blocks.
!!
!!***

subroutine Grid_receiveInputData(localNumBlocks, alnblocks, xx)
  implicit none
  integer, intent(IN) :: localNumBlocks, alnblocks, xx
end subroutine Grid_receiveInputData
