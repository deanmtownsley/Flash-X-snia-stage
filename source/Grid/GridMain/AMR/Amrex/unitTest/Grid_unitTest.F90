!!****if* source/Grid/unitTest/Grid_unitTest
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
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  call Grid_unitTest(integer(in):: fileUnit,
!!                     logical(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This unit test exercises the data accessing functions of the Grid unit
!!  The routine has direct access to all the mesh data structures such as 
!!  "unk", "facevarx" etc. It uses the Grid_getBlk/Point/RowData functions 
!!  to fetch some or all of the block data, and then compares it with
!!  the corresponding section of the appropriate array.
!!
!! ARGUMENTS
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!!
!!***

subroutine Grid_unitTest(fileUnit,perfect)
  implicit none
  
  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors

  return
end subroutine Grid_unitTest
