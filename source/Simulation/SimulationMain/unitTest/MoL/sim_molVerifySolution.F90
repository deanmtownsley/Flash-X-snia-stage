!!****f* source/Simulation/SimulationMain/unitTest/MoL/sim_molVerifySolution
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
!!  sim_molVerifySolution
!!
!!
!! SYNOPSIS
!!  call sim_molVerifySolution(real,    intent(in)  :: t,
!!                             real,    intent(in)  :: dt,
!!                             logical, intent(out) :: valid,
!!                             real,    intent(out) :: maxError)
!!
!!
!! DESCRIPTION
!!    Verify the solution at the specified time.  An implementation
!!    may utilize an analytic/baseline solution or some other
!!    error-estimator to determine if the current solution is valid
!!    for the purposes of passing a unit test.
!!
!! ARGUMENTS
!!
!!    t        : The current time that the solution is at
!!    dt       : The current timestep
!!    valid    : Is this a valid solution
!!    maxError : The maximum error present in the solution that was utilized
!!               to determine if the solution was valid
!!
!!***
subroutine sim_molVerifySolution(t, dt, valid, maxError)
   implicit none

   real, intent(in) :: t, dt
   logical, intent(out) :: valid
   real, intent(out) :: maxError

   ! The stub should break the unit test...
   valid = .false.
   maxError = 1d40
end subroutine sim_molVerifySolution
