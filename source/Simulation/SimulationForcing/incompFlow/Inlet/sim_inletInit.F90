!!***if* source/Simulation/SimulationForcing/incompFlow/Inlet/sim_inletInit
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

subroutine sim_inletInit()

   use sim_inletData
   use Grid_interface, ONLY: Grid_getDomainBC
   use Simulation_data, ONLY: sim_meshMe
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use Driver_interface, ONLY: Driver_abort

   implicit none
   integer :: idimn, ibound, domainBC(LOW:HIGH, MDIM)

   call RuntimeParameters_get('sim_inletSink', sim_inletSink)
   call RuntimeParameters_get('sim_inletBuffer', sim_inletBuffer)
   call RuntimeParameters_get('sim_inletGrowthRate', sim_inletGrowthRate)
   call Grid_getDomainBC(domainBC)

   sim_inletFlag = 0

   do idimn = 1, NDIM
      do ibound = LOW, HIGH

         select case (domainBC(ibound, idimn))

         case (INFLOW_INS)
            sim_inletFlag(ibound, idimn) = 1

         end select
      end do
   end do

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'sim_inletSink = ', sim_inletSink
      write (*, *) 'sim_inletBuffer =', sim_inletBuffer
      write (*, *) 'sim_inletGrowthRate =', sim_inletGrowthRate
      write (*, *) 'Inlet Flag Low  =', sim_inletFlag(LOW, :)
      write (*, *) 'Inlet Flag High =', sim_inletFlag(HIGH, :)
   end if

end subroutine sim_inletInit
