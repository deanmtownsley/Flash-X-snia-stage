!!***if* source/Simulation/SimulationMain/incompFlow/PoolBoiling/localAPI/heater/sim_heaterTagSites
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

subroutine sim_heaterTagSites(stime)

   use Simulation_data, ONLY: sim_meshMe
   use sim_heaterData

   implicit none
   include "Flashx_mpi.h"
   real, intent(in) :: stime

   integer :: htr,ierr,isite
   type(sim_heaterType), pointer :: heater
   integer :: TA(2),count_rate
   real*8  :: ET

   CALL SYSTEM_CLOCK(TA(1),count_rate)
 
   do htr=1,sim_numHeaters

    heater => sim_heaterInfo(htr)

    do isite=1,heater%numSites

       call MPI_Allreduce(heater%siteIsAttachedCurr(isite), heater%siteIsAttachedCurr(isite), &
                          1, FLASH_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

       if(heater%siteIsAttachedPrev(isite) .eqv. .true.  .and. &
          heater%siteIsAttachedCurr(isite) .eqv. .false.) heater%siteTimeStamp(isite) = stime

       if (sim_meshMe .eq. MASTER_PE .and. sim_heaterShowInfo) &
           write(*,'(A,I2,A,I3,A,L1,A,2g14.6)')&
                   ' Heater:',htr,&
                   ' Site:',isite,&
                   ' IsAttached:',heater%siteIsAttachedCurr(isite),&
                   ' TimeStamp:',heater%siteTimeStamp(isite)

    end do

    heater%siteIsAttachedPrev = heater%siteIsAttachedCurr
    heater%siteIsAttachedCurr = .false.

   end do

   CALL SYSTEM_CLOCK(TA(2),count_rate)
   ET=REAL(TA(2)-TA(1))/count_rate
   if (sim_meshMe .eq. MASTER_PE)  write(*,*) 'Total sim_heater TagSites Time =',ET

   return

end subroutine sim_heaterTagSites
