!!***if* source/physics/sourceTerms/Heater/HeaterMain/Heater_tagSites
!!
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
!!
!!***
#include "constants.h"
#include "Simulation.h"

subroutine Heater_tagSites(stime)

   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Heater_type, ONLY: Heater_type_t
   use Heater_data

   implicit none
   include "Flashx_mpi.h"
   real, intent(in) :: stime

   integer :: htr, ierr, isite
   type(Heater_type_t), pointer :: heater

   call Timers_start("Heater_tagSites")

#ifdef MULTIPHASE_EVAPORATION
   do htr = 1, htr_numHeaters

      heater => htr_heaterInfo(htr)

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

         if (htr_meshMe .eq. MASTER_PE .and. htr_showInfo) &
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

   call Timers_stop("Heater_tagSites")

   return

end subroutine Heater_tagSites
