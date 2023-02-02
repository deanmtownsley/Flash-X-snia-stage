!!****if* source/Simulation/SimulationMain/UnshiftedGaugeWave/Simulation_init
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
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  call Simulation_init()
!!
!!
!! DESCRIPTION
!!  Initializes all the parameters needed for a particular simulation
!!
!! ARGUMENTS
!!
!!
!!***
subroutine Simulation_init()
   use Simulation_data, only: sim_dt

   use RuntimeParameters_interface, only: RuntimeParameters_get

   implicit none

   call RuntimeParameters_get("sim_dt", sim_dt)

end subroutine Simulation_init
