!!****f* source/Grid/Grid_unitTest
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
!!  Grid_unitTest(integer, intent(in):: fileUnit,
!!                logical, intent(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  The Grid unit test has several implementations. There is a general test,
!!  which exercises the data access interface of the Grid unit, for example
!!  the get/putData functions, and can be used with all GridMain
!!  implementations. There is also a test specific to
!!  the Uniform Grid which tests internal working of the 
!!  UG such as parallel data exchanges, domain setup, etc. 
!!   
!!
!! ARGUMENTS
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!! NOTES
!!
!!***

subroutine Grid_unitTest(fileUnit,perfect)

  implicit none

  
  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors

  perfect = .false.

  return
 
end subroutine Grid_unitTest
