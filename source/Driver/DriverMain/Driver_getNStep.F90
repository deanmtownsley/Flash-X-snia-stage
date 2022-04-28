!!****if* source/Driver/DriverMain/Driver_getNStep
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
!!  Driver_getNStep
!!
!! SYNOPSIS
!!  
!!
!!  Driver_getNStep(real(out) :: nstep)
!!  
!! DESCRIPTION 
!!
!!  Accessor function that returns the current step number in the
!! simulation from the Driver. 
!!
!! ARGUMENTS
!!  nstep - returned value, current step number 
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
!!***

subroutine Driver_getNStep(nstep)

  use Driver_data, ONLY : dr_nStep

implicit none
  integer, intent(out) :: nstep
  
  nstep = dr_nStep

end subroutine Driver_getNStep

