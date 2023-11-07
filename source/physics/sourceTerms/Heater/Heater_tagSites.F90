!!***if* source/Simulation/SimulationForcing/incompFlow/Heater/Heater_TagSites
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

subroutine Heater_TagSites(stime)

   use Simulation_data, ONLY: sim_meshMe
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Heater_Data

   implicit none
   include "Flashx_mpi.h"
   real, intent(in) :: stime

   integer :: htr, ierr, isite
   type(Heater_Type), pointer :: heater

   call Timers_start("Heater_TagSites")

#ifdef MULTIPHASE_EVAPORATION
   do htr = 1, sim_numHeaters

      heater => Heater_Info(htr)

      !! DEVNOTE (10/20/2023): This all reduce is not relevant anymore since implementation
      !!                       of Heater_MapSitesToProc. Leaving it here for legacy
      !!
      !call Timers_start("consolidate site status")
      !call MPI_Allreduce(MPI_IN_PLACE, heater%siteIsAttachedCurr, &
      !                   heater%numSites, FLASH_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
      !call Timers_stop("consolidate site status")

      do isite = 1, heater%numSitesProc

         if (heater%siteIsAttachedPrev(isite) .eqv. .true. .and. &
             heater%siteIsAttachedCurr(isite) .eqv. .false.) heater%siteTimeStamp(isite) = stime

         if (sim_meshMe .eq. MASTER_PE .and. Heater_ShowInfo) &
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

   call Timers_stop("Heater_TagSites")

   return

end subroutine Heater_TagSites
