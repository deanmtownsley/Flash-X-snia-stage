!!****f* source/Simulation/SimulationMain/incompFlow/FlowBoiling/Simulation_finalize
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
!! NAME
!!  Simulation_finalize
!!
!! SYNOPSIS
!!
!!  Simulation_finalize()
!!
!! DESCRIPTION
!!
!!  This dummy function cleans up the Simulation unit, deallocates memory, etc.
!!  However, as nothing needs to be done, only this stub is included.
!!
!! ARGUMENTS
!!
!!
!!
!!***

subroutine Simulation_finalize()

   use sim_heaterInterface, ONLY: sim_heaterFinalize
   use sim_outletInterface, ONLY: sim_outletFinalize
   use sim_inletInterface, ONLY: sim_inletFinalize

   implicit none

   call sim_heaterFinalize()
   call sim_inletFinalize()
   call sim_outletFinalize()

end subroutine Simulation_finalize
