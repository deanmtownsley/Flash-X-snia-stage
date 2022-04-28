!!****f* source/Multispecies/Multispecies_unitTest
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
!!  Multispecies_unitTest
!!
!! SYNOPSIS
!!
!!  Multispecies_unitTest(integer, intent(in)::fileUnit,
!!                        logical, intent(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This is a unitTest setup for testing the Multispecies
!!      unit.  Normally called from Driver_evolveAll within a Simulation unitTest
!!  See, for example, source/Simulation/SimulationMain/unitTest/Multispecies
!!
!! ARGUMENTS
!!
!!    fileUnit  -- number of file unit for diagnostic output
!!    perfect   -- if .true., unitTest has returned correctly
!! 
!! NOTES
!!  In the Simulation unit, you must set your Config file to 
!!    REQUIRES Multispecies/MultispeciesMain/unitTest
!!
!!***

subroutine Multispecies_unitTest(fileUnit,perfect)
!------------------------------------------------------------------------
    implicit none
    
 
    integer, intent(in)         :: fileUnit
    logical, intent(inout)      :: perfect
 
end subroutine Multispecies_unitTest
