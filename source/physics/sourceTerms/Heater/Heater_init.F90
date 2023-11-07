!!***if* source/Simulation/SimulationForcing/incompFlow/Heater/Heater_Init
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

subroutine Heater_Init()

   use Simulation_data, ONLY: sim_meshMe
   use Heater_Data
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use Heater_Interface, ONLY: Heater_Read

   implicit none
   character(len=30) :: heaterFile
   integer           :: htr

   call RuntimeParameters_get('sim_numHeaters', sim_numHeaters)
   call RuntimeParameters_get('sim_nucSeedRadius', sim_nucSeedRadius)
   call RuntimeParameters_get('Heater_Name', Heater_Name)
   call RuntimeParameters_get('Heater_ShowInfo', Heater_ShowInfo)
#ifdef SIM_HEATER_ANN_SEARCH
   call RuntimeParameters_get("Heater_AnnQueries", Heater_AnnQueries)
#endif

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'sim_numHeaters=', sim_numHeaters
      write (*, *) 'sim_nucSeedRadius=', sim_nucSeedRadius
#ifdef SIM_HEATER_ANN_SEARCH
      write (*, *) 'Heater_AnnQueries=', Heater_AnnQueries
#endif
   end if

#ifdef SIM_HEATER_ANN_SEARCH
   allocate (Heater_AnnIdx(Heater_AnnQueries))
#endif

   allocate (Heater_Info(sim_numHeaters))

   do htr = 1, sim_numHeaters
      write (heaterFile, "(A,A,I4.4)") trim(Heater_Name), '_hdf5_htr_', htr
      call Heater_Read(htr, heaterFile)
   end do

end subroutine Heater_Init
