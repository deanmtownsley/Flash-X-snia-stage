!!****if* source/Grid/localAPI/gr_setBlockType
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
!!  gr_setBlockType
!!
!! SYNOPSIS
!!
!!  gr_setBlockType(integer(IN)  :: blockID,
!!                     integer(IN)  :: type)
!!  
!! DESCRIPTION 
!!  
!!
!! ARGUMENTS 
!!
!!  blockID : the local blockid
!!  type    : the type to which block is to be set
!!
!!
!!
!! NOTES
!!
!!  gr_setBlockType is a mutator function that sets a block to the specified Type
!!  For example when tracing the path of a ray, a blocktype can be set to "TRAVERSE"
!!  if the ray passed through it.
!!
!!
!!***

subroutine gr_setBlockType(blockID,type)
  implicit none

  integer, intent(IN) :: blockID, type

end subroutine gr_setBlockType
