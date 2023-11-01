!!***if* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterFinalize
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

subroutine sim_heaterFinalize()

   use sim_heaterData, ONLY: sim_heaterInfo
#ifdef SIM_HEATER_ANN_SEARCH
   use sim_heateData, ONLY: sim_heaterAnnIdx
#endif

   implicit none

   deallocate (sim_heaterInfo)
#ifdef SIM_HEATER_ANN_SEARCH
   deallocate (sim_heaterAnnIdx)
#endif

end subroutine sim_heaterFinalize
