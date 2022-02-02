!!****if* source/Driver/DriverMain/Driver_getSimTime
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
!!  Driver_getSimTime
!!
!! SYNOPSIS
!!  
!!
!!  call Driver_getSimTime(real(out) :: simulationTime,
!!                OPTIONAL,integer(out) :: simGeneration)
!!  
!! DESCRIPTION 
!!
!!  Accessor function that returns the current simulation time from the Driver
!!
!! ARGUMENTS
!!  simulationTime - returned value, current run simulation time
!!  simGeneration  - returned value, current generation of the solution
!!                   at this simulation time
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!
!!***

subroutine Driver_getSimTime(simulationTime, simGeneration)

  use Driver_data, ONLY : dr_simTime, dr_simGeneration

  implicit none
  real, intent(out) :: simulationTime
  integer, intent(out), OPTIONAL :: simGeneration
  
  simulationTime = dr_simTime
  if (present(simGeneration)) simGeneration = dr_simGeneration

end subroutine Driver_getSimTime

