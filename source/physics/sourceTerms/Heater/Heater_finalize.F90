!!***if* source/Simulation/SimulationForcing/incompFlow/Heater/Heater_Finalize
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

#include "Simulation.h"

subroutine Heater_Finalize()

   use Heater_Data, ONLY: Heater_Info
#ifdef SIM_HEATER_ANN_SEARCH
   use sim_heateData, ONLY: Heater_AnnIdx
#endif

   implicit none

   deallocate (Heater_Info)
#ifdef SIM_HEATER_ANN_SEARCH
   deallocate (Heater_AnnIdx)
#endif

end subroutine Heater_Finalize
