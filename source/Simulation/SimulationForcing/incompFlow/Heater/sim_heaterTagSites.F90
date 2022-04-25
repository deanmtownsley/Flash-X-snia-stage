!!***if* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterTagSites
!!
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
!!
!!***
#include "constants.h"
#include "Simulation.h"

subroutine sim_heaterTagSites(stime)

   use Simulation_data, ONLY: sim_meshMe
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use sim_heaterData

   implicit none
   include "Flashx_mpi.h"
   real, intent(in) :: stime

   integer :: htr, ierr, isite
   type(sim_heaterType), pointer :: heater

   call Timers_start("sim_heaterTagSites")

#ifdef MULTIPHASE_EVAPORATION
   do htr = 1, sim_numHeaters

      heater => sim_heaterInfo(htr)

      call Timers_start("consolidate site status")
      call MPI_Allreduce(MPI_IN_PLACE, heater%siteIsAttachedCurr, &
                         heater%numSites, FLASH_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
      call Timers_stop("consolidate site status")

      do isite = 1, heater%numSites

         if (heater%siteIsAttachedPrev(isite) .eqv. .true. .and. &
             heater%siteIsAttachedCurr(isite) .eqv. .false.) heater%siteTimeStamp(isite) = stime

         if (sim_meshMe .eq. MASTER_PE .and. sim_heaterShowInfo) &
            write (*, '(A,I2,A,I3,A,L1,A,2g14.6)') &
            ' Heater:', htr, &
            ' Site:', isite, &
            ' IsAttached:', heater%siteIsAttachedCurr(isite), &
            ' TimeStamp:', heater%siteTimeStamp(isite)

      end do

      heater%siteIsAttachedPrev = heater%siteIsAttachedCurr
      heater%siteIsAttachedCurr = .false.

   end do
#endif

   call Timers_stop("sim_heaterTagSites")

   return

end subroutine sim_heaterTagSites
