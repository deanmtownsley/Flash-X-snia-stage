!!****if* source/Simulation/SimulationMain/Sedov/Simulation_adjustEvolution
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
!!  Simulation_adjustEvolution( integer(IN) :: blkcnt,
!!                              integer(IN) :: blklst(blkcnt),
!!                              integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!  
!! ARGUMENTS
!!  blkcnt - number of blocks
!!  blklst - block list
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***

#include "Simulation.h"

subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)
  use Simulation_interface, ONLY: Simulation_computeAnalytical
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real, intent(in) :: dt
  real, intent(in) :: stime

#if NDIM==1
  real    :: tnew
  integer :: lb

  call Timers_start("adjustEvo")


!!$  ! We want the NEW time that we are advancing TO.
!!$  ! This is based on where Simulation_adjustEvolution is called from Driver_evolveAll:
!!$  ! dr_simTime has not been updated yet at this point.
!!$  !!DEV: Check whether this needs adjustments for non-"Unsplit" Driver_evolveAll.F90 !
!!$  !!DEV: Check whether this needs adjustments for when STS is used!
!!$  tnew = stime+dt
!!$
!!$  do lb = 1, blkcnt
!!$     call Simulation_computeAnalytical(blklst(lb), tnew)
!!$  end do
  call Timers_stop("adjustEvo")
#endif

end subroutine Simulation_adjustEvolution
