!!****if* source/Simulation/SimulationMain/Brusselator/Driver_computeDt
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  Driver_computeDt
!!
!!
!! SYNOPSIS
!!
!!  Driver_computeDt(integer(IN) :: nbegin,
!!                  integer(IN) :: nstep,
!!                  real(IN)    :: simTime,
!!                  real(IN)    :: dtOld,
!!                  real(OUT)   :: dtNew)
!!
!! DESCRIPTION
!!
!!  Determine the stability-limited time step.
!!  This timestep is determined using information from the included
!!  physics modules - many different timestep limiters are polled.
!!
!!  The global driver might use a different (hopefully smaller) time
!!  step, to match a file write time (tplot or trstr) or if the
!!  simulation end time has been reached; such possibilities are
!!  not considered here.
!!
!! ARGUMENTS
!!  nbegin - first step of the simulation (this is only used
!!              to determine if a label header should be written to
!!              the screen)
!!  nstep - current step of the simulation
!!  simTime - current simulation time of the run
!!  dtOld - the dt from the timestep that we just finished
!!         (it's old because we be using dtOld to calculate
!!          and return the dt for the next timestep (dtNew)
!!  dtNew - returned value of the dt calculated for the next timestep
!!
!!
!!***
subroutine Driver_computeDt(nbegin, nstep, simTime, dtOld, dtNew)
   use Simulation_data, only: sim_k

   implicit none

   integer, intent(in)  :: nbegin, nstep
   real, intent(in)  :: simTime    !! current simulation time
   real, intent(in)  :: dtOld      !! last time step we used
   real, intent(out) :: dtNew      !! the new timestep we get. to be returned.

   dtNew = 0.1/2.0**sim_k
end subroutine Driver_computeDt
