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
#include "Simulation.h"

subroutine sim_outletInit()

   use sim_outletData
   use Grid_interface, ONLY: Grid_getDomainBC
   use Simulation_data, ONLY: sim_meshMe
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use Driver_interface, ONLY: Driver_abort

   integer :: idimn, ibound, domainBC(LOW:HIGH, MDIM)
   character(len=100) :: errorMessage

   call RuntimeParameters_get('sim_outletSink', sim_outletSink)
   call RuntimeParameters_get('sim_outletBuffer', sim_outletBuffer)
   call RuntimeParameters_get('sim_outletGrowthRate', sim_outletGrowthRate)
   call Grid_getDomainBC(domainBC)

   sim_outletFlag = 0
   sim_QOut = 0.
   sim_QOutLiq = 0.
   sim_QOutGas = 0.

   do idimn = 1, NDIM
      do ibound = LOW, HIGH
         select case (domainBC(ibound, idimn))
         case (OUTFLOW_INS)
            sim_outletFlag(ibound, idimn) = 1
         end select
      end do
   end do

   if (sim_meshMe .eq. MASTER_PE) then
      write (*, *) 'sim_outletSink=', sim_outletSink
      write (*, *) 'sim_outletBuffer=', sim_outletBuffer
      write (*, *) 'sim_outletGrowthRate=', sim_outletGrowthRate
      write (*, *) 'Outlet Flag Low  =', sim_outletFlag(LOW, :)
      write (*, *) 'Outlet Flag High =', sim_outletFlag(HIGH, :)
   end if

end subroutine sim_outletInit
