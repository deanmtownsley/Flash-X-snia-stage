!!***if* source/physics/sourceTerms/Inlet/InletMain/Inlet_init
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

subroutine Inlet_init()

   use Inlet_data
   use Grid_interface, ONLY: Grid_getDomainBC
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use Driver_interface, ONLY: Driver_abort, Driver_getMyPe

   implicit none
   integer :: idimn, ibound, domainBC(LOW:HIGH, MDIM)

   call Driver_getMype(MESH_COMM, in_meshMe)

   call RuntimeParameters_get('in_sink', in_sink)
   call RuntimeParameters_get('in_buffer', in_buffer)
   call RuntimeParameters_get('in_growthRate', in_growthRate)
   call Grid_getDomainBC(domainBC)

   in_flag = 0

   do idimn = 1, NDIM
      do ibound = LOW, HIGH

         select case (domainBC(ibound, idimn))

         case (INFLOW_INS)
            in_flag(ibound, idimn) = 1

         end select
      end do
   end do

   if (in_meshMe .eq. MASTER_PE) then
      write (*, *) 'in_sink = ', in_sink
      write (*, *) 'in_buffer =', in_buffer
      write (*, *) 'in_growthRate =', in_growthRate
      write (*, *) 'Inlet Flag Low  =', in_flag(LOW, :)
      write (*, *) 'Inlet Flag High =', in_flag(HIGH, :)
   end if

end subroutine Inlet_init
