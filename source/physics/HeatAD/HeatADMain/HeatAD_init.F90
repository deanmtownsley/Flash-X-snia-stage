subroutine HeatAD_init(restart)
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

   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use Driver_interface, ONLY: Driver_getMype, Driver_getNumProcs, Driver_getComm
   use ht_advectionInterface, ONLY: ht_advectionInit
   use HeatAD_data

#include "constants.h"
#include "Simulation.h"
#include "HeatAD.h"

   implicit none
   logical, INTENT(IN) :: restart
   logical :: useIncompNS, useMultiphase

   call RuntimeParameters_get("useHeatAD", ht_useHeatAD)
   if (.NOT. ht_useHeatAD) RETURN

   call Driver_getMype(MESH_COMM, ht_meshMe)
   call Driver_getNumProcs(MESH_COMM, ht_meshNumProcs)
   call Driver_getComm(MESH_COMM, ht_meshComm)

   call RuntimeParameters_get('ht_Twall_high', ht_Twall_high)
   call RuntimeParameters_get('ht_Twall_low', ht_Twall_low)
   call RuntimeParameters_get('ht_Prandtl', ht_Prandtl)
   call RuntimeParameters_get('ht_Tbulk', ht_Tbulk)

   ht_invReynolds = 1.

   call RuntimeParameters_get('useIncompNS', useIncompNS)
   if (useIncompNS) call RuntimeParameters_get('ins_invReynolds', ht_invReynolds)

   ht_thcoGas = 1.
   ht_CpGas = 1.

   call RuntimeParameters_get('useMultiphase', useMultiphase)
   if (useMultiphase) then
      call RuntimeParameters_get('mph_thcoGas', ht_thcoGas)
      call RuntimeParameters_get('mph_cpGas', ht_CpGas)
      call RuntimeParameters_get('mph_Tsat', ht_Tsat)
   end if

   if (ht_meshMe .eq. MASTER_PE) then
      write (*, *) 'ht_Twall_high   =', ht_Twall_high
      write (*, *) 'ht_Twall_low    =', ht_Twall_low
      write (*, *) 'ht_Prandtl      =', ht_Prandtl
      write (*, *) 'ht_Tsat         =', ht_Tsat
      write (*, *) 'ht_Tbulk        =', ht_Tbulk
      write (*, *) 'ht_invReynolds  =', ht_invReynolds
      write (*, *) 'ht_thcoGas      =', ht_thcoGas
      write (*, *) 'ht_CpGas        =', ht_CpGas
   end if

   call ht_advectionInit()

end subroutine HeatAD_init
