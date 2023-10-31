!!***if* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterInit
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

subroutine sim_heaterInit()

   use Simulation_data, ONLY: sim_meshMe
   use sim_heaterData
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use sim_heaterInterface, ONLY: sim_heaterRead

   implicit none
   character(len=30) :: heaterFile
   integer           :: htr

   call RuntimeParameters_get('sim_numHeaters', sim_numHeaters)
   call RuntimeParameters_get('sim_nucSeedRadius', sim_nucSeedRadius)
   call RuntimeParameters_get('sim_heaterName', sim_heaterName)
   call RuntimeParameters_get('sim_heaterShowInfo', sim_heaterShowInfo)
#ifdef SIM_HEATER_ANN_SEARCH
   call RuntimeParameters_get("sim_heaterAnnQueries", sim_heaterAnnQueries)
#endif

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'sim_numHeaters=', sim_numHeaters
      write (*, *) 'sim_nucSeedRadius=', sim_nucSeedRadius
#ifdef SIM_HEATER_ANN_SEARCH
      write (*, *) 'sim_heaterAnnQueries=', sim_heaterAnnQueries
#endif
   end if

#ifdef SIM_HEATER_ANN_SEARCH
   allocate (sim_heaterAnnIdx(sim_heaterAnnQueries))
#endif

   allocate (sim_heaterInfo(sim_numHeaters))

   do htr = 1, sim_numHeaters
      write (heaterFile, "(A,A,I4.4)") trim(sim_heaterName), '_hdf5_htr_', htr
      call sim_heaterRead(htr, heaterFile)
   end do

end subroutine sim_heaterInit
