!!****f* source/physics/Gravity/Gravity_unitTest
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
!! NAME
!!
!!  Gravity_unitTest
!! 
!! SYNOPSIS
!!
!!  call Gravity_unitTest()
!!                    integer(IN) :: fileUnit,
!!                    logical(OUT) :: perfect
!!
!! DESCRIPTION
!!
!! This function is the unit test for the Gravity unit. It is invoked in
!! the setup unitTest/Gravity. The Config file for Gravity unit test setup
!! requests an extra variable in the main grid data structure for
!! storing the analytical solution
!!
!!  ARGUMENTS 
!!   
!!   
!! 
!!   fileUnit : unit number for file opened by the unitTest/Gravity setup
!!              in which to write results of the test
!!
!!   perfect : indicates test ran without error is true.
!!
!!  PARAMETERS
!!
!!  eintSwitch  a rarely used switch which ensures that internal energy calculations 
!!        maintain sufficient precision. Important only if energyTotal is dominated 
!!        by energyKinetic.
!!
!!***

subroutine Gravity_unitTest( fileUnit, perfect)

  use Logfile_interface, ONLY:  Logfile_stampMessage

  implicit none

  integer, intent(in) :: fileUnit
  logical, intent(out) :: perfect

  perfect = .false.

  write(*,*)"[Gravity_unitTest] Failure.  This is only a stub."
  call Logfile_stampMessage("[Gravity_unitTest] Failure.  This is only a stub.")

  return
end subroutine Gravity_unitTest




