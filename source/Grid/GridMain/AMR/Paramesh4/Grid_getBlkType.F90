!!****if* source/Grid/GridMain/paramesh/Grid_getBlkType
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
!!  Grid_getBlkType
!!
!! SYNOPSIS
!!
!!
!!  Grid_getBlkType(integer(IN)  :: blockID,
!!                  integer(OUT) :: blkType)
!!  
!! DESCRIPTION 
!!  Get the type of block. 
!!  Valid types are defined in constants.h, and include
!!  LEAF, PARENT and ANCESTOR for paramesh.
!!
!! ARGUMENTS
!!  blockID - the local block number
!!  blkType - returned value
!!
!!***

subroutine Grid_getBlkType(blockID, blkType)
  use tree, ONLY : nodetype
  implicit none
#include "constants.h"
  integer,intent(in) :: blockID
  integer,intent(out) :: blkType
  blkType=nodetype(blockID)
  return
end subroutine Grid_getBlkType














