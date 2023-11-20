!!***if* source/physics/sourceTerms/Outlet/OutletMain/Outlet_init
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

subroutine Outlet_init()

   use Outlet_data
   use Grid_interface, ONLY: Grid_getDomainBC
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use Driver_interface, ONLY: Driver_abort

   integer :: idimn, ibound, domainBC(LOW:HIGH, MDIM)
   character(len=100) :: errorMessage

   call Driver_getMype(MESH_COMM, out_meshMe)

   call RuntimeParameters_get('xmin', out_xMin)
   call RuntimeParameters_get('ymin', out_yMin)
   call RuntimeParameters_get('xmax', out_xMax)
   call RuntimeParameters_get('ymax', out_yMax)
   call RuntimeParameters_get('zmin', out_zmin)
   call RuntimeParameters_get('zmax', out_zmax)

   call RuntimeParameters_get('out_sink', out_sink)
   call RuntimeParameters_get('out_buffer', out_buffer)
   call RuntimeParameters_get('out_growthRate', out_growthRate)
   call RuntimeParameters_get('out_velRefScale', out_velRefScale)

   call Grid_getDomainBC(domainBC)

   out_flag = 0
   out_QOut = 0.
   out_QOutLiq = 0.
   out_QOutGas = 0.

   do idimn = 1, NDIM
      do ibound = LOW, HIGH
         select case (domainBC(ibound, idimn))
         case (OUTFLOW_INS, EXTRAP_INS)
            out_flag(ibound, idimn) = 1
         end select
      end do
   end do

   if (out_meshMe .eq. MASTER_PE) then
      write (*, *) 'out_sink=', out_sink
      write (*, *) 'out_buffer=', out_buffer
      write (*, *) 'out_growthRate=', out_growthRate
      write (*, *) 'out_velRefScale=', out_velRefScale
      write (*, *) 'Outlet Flag Low  =', out_flag(LOW, :)
      write (*, *) 'Outlet Flag High =', out_flag(HIGH, :)
   end if

end subroutine Outlet_init
