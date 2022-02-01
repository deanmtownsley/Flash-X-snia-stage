!!****f* source/Simulation/Simulation_adjustEvolution
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
!!  Simulation_adjustEvolution
!!
!!
!! SYNOPSIS
!!  Simulation_adjustEvolution( integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!  
!! ARGUMENTS
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***
subroutine Simulation_adjustEvolution(nstep, dt, stime)
  implicit none

  integer, intent(in) :: nstep
  real, intent(in) :: dt
  real, intent(in) :: stime

end subroutine Simulation_adjustEvolution
