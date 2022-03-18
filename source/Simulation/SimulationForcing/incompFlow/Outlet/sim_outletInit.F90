!!***if* source/Simulation/SimulationForcing/incompFlow/Outlet/sim_outletInit
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

subroutine sim_outletInit()

   use Grid_interface, ONLY: Grid_getDomainBC
   use Simulation_data, ONLY: sim_meshMe
   use sim_outletData
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get

   implicit none

   call Grid_getDomainBC(sim_domainBC)

   call RuntimeParameters_get('sim_outletSink', sim_outletSink)
   call RuntimeParameters_get('sim_outletBuffer', sim_outletBuffer)
   call RuntimeParameters_get('sim_outletGrowthRate', sim_outletGrowthRate)

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'sim_outletSink=', sim_outletSink
      write (*, *) 'sim_outletBuffer=', sim_outletBuffer
      write (*, *) 'sim_outletGrowthRate=', sim_outletGrowthRate
   end if

   sim_outletVel = 0.

end subroutine sim_outletInit
