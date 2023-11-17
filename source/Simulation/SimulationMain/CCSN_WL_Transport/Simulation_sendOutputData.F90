!!****if* source/Simulation/SimulationMain/CCSN_WL_Transport/Simulation_sendOutputData
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
!!  Simulation_sendOutputData
!!
!! SYNOPSIS
!!
!!  Simulation_sendOutputData()
!!
!! DESCRIPTION
!!
!! This routine sends the scalar variables owned by the Simulation unit
!! to the IO unit, to be written to a checkpoint file.
!!
!!
!!***
subroutine Simulation_sendOutputData()
   use Simulation_data, only: sim_postBounce, sim_bounceTime
   use IO_interface, only: IO_setScalar
   use Driver_interface, only: Driver_getSimTime
   ! use Eos_wlInterface, only: Eos_wlDetectBounce

   implicit none

   real :: postBounceTime, simTime

   if (.not. sim_postBounce) &
      call sim_detectBounce(sim_postBounce, sim_bounceTime)

   call IO_setScalar("postBounce", sim_postBounce)

   postBounceTime = 0.0
   if (sim_postBounce) then
      call IO_setScalar("bounceTime", sim_bounceTime)

      call Driver_getSimTime(simTime)
      postBounceTime = simTime - sim_bounceTime
   end if

   call IO_setScalar("postBounceTime", postBounceTime)

end subroutine Simulation_sendOutputData

