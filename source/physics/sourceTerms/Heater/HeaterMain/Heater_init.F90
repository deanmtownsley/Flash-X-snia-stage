!!***if* source/physics/sourceTerms/Heater/HeaterMain/Heater_init
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

subroutine Heater_init()

   use Heater_data
   use Driver_interface, ONLY: Driver_getMype
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use htr_interface, ONLY: htr_read

   implicit none
   character(len=30) :: heaterFile
   integer           :: heater

   call Driver_getMype(MESH_COMM, htr_meshMe)

   call RuntimeParameters_get('htr_numHeaters', htr_numHeaters)
   call RuntimeParameters_get('htr_nucSeedRadius', htr_nucSeedRadius)
   call RuntimeParameters_get('htr_heaterName', htr_heaterName)
   call RuntimeParameters_get('htr_showInfo', htr_showInfo)
#ifdef HEATER_ANN_SEARCH
   call RuntimeParameters_get("htr_annQueries", htr_annQueries)
#endif

   if (htr_meshMe .eq. MASTER_PE) then
      write (*, *) 'htr_numHeaters=', htr_numHeaters
      write (*, *) 'htr_nucSeedRadius=', htr_nucSeedRadius
#ifdef HEATER_ANN_SEARCH
      write (*, *) 'htr_annQueries=', htr_annQueries
#endif
   end if

#ifdef HEATER_ANN_SEARCH
   allocate (htr_annIdx(htr_annQueries))
#endif

   allocate (htr_heaterInfo(htr_numHeaters))

   do heater = 1, htr_numHeaters
      write (heaterFile, "(A,A,I4.4)") trim(htr_heaterName), '_hdf5_htr_', heater
      call htr_read(heater, heaterFile)
   end do

end subroutine Heater_init
