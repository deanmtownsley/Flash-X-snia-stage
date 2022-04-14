!!****if* source/physics/IncompNS/IncompNSMain/ins_init
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
!!
!!*****
!!
#include "constants.h"
#include "Simulation.h"

subroutine ins_init()

   use IncompNS_data
   use Driver_interface, ONLY: Driver_abort
   use Grid_interface, ONLY: GRID_PDE_BND_NEUMANN, GRID_PDE_BND_DIRICHLET, &
                             Grid_getDomainBoundBox, Grid_getDomainBC

   implicit none
   integer :: idimn, ibound, eachBoundary

   call Grid_getDomainBoundBox(ins_globalDomain)
   call Grid_getDomainBC(ins_domainBC)

   ins_pressureBC_types(:) = GRID_PDE_BND_NEUMANN
   ins_pressureBC_values(:, :) = 0.0

   do idimn = 1, NDIM
      do ibound = LOW, HIGH

         eachBoundary = 2*(idimn - 1) + ibound

         select case (ins_domainBC(ibound, idimn))
         case (PERIODIC)
            ins_pressureBC_types(eachBoundary) = GRID_PDE_BND_NEUMANN

         case (SLIP_INS, NOSLIP_INS, INFLOW_INS, MOVLID_INS, EXTRAP_INS)
            ins_pressureBC_types(eachBoundary) = GRID_PDE_BND_NEUMANN

         case (OUTFLOW_INS)
            ins_pressureBC_types(eachBoundary) = GRID_PDE_BND_DIRICHLET

         case default
            if (ins_meshMe .eq. MASTER_PE) then
               write (*, *) 'IncompNS Error: Boundary Conditions match for Poisson Solver not defined.'
               write (*, *) 'IncompNS Error: LOW-HIGH,AXIS=', ibound, idimn
               write (*, *) 'IncompNS Error: ins_domainBC(ibound,idimn) =', ins_domainBC(ibound, idimn)
            end if
            call Driver_abort('IncompNS Error: BCs do not have matching Poisson solver BCs')
         end select
      end do
   end do

   ins_poisfact = 1.0

   if (ins_meshMe .eq. MASTER_PE) write (*, *) 'IncompNS: Pressure Boundary Conditions'
   do idimn = 1, NDIM
   do ibound = LOW, HIGH
      eachBoundary = 2*(idimn - 1) + ibound
      if (ins_meshMe .eq. MASTER_PE) then
         write (*, *) 'LOW-HIGH,AXIS=', ibound, idimn
         write (*, *) 'ins_domainBC(ibound,idimn), ins_pressureBC(ibound,idimn) =', ins_domainBC(ibound, idimn), &
            ins_pressureBC_types(eachBoundary)
      end if
   end do
   end do

end subroutine ins_init
