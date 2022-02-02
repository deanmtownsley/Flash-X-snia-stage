!!****f* source/PhysicalConstants/PhysicalConstants_unitTest
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
!!  PhysicalConstants_unitTest
!!
!! SYNOPSIS
!!
!!  PhysicalConstants_unitTest(integer, intent(in)::fileUnit,
!!                             logical, intent(out)::perfect  )
!!
!! DESCRIPTION
!!
!!  This is a unitTest setup for testing the PhysicalConstants
!!      unit.  Normally called from Driver_evolveAll within a Simulation unitTest
!!  See, for example, source/Simulation/SimulationMain/unitTest/PhysConst
!!
!! ARGUMENTS
!!      
!!     fileUnit - integer(in)   file for diagnostic output
!!     perfect  - logical(out)  TRUE if all tests are passed
!!
!! NOTES
!!
!!***

subroutine PhysicalConstants_unitTest(fileUnit,perfect)

  implicit none

  integer, intent(in)         :: fileUnit
  logical, intent(out)        :: perfect
  
  return 
  
end subroutine PhysicalConstants_unitTest
