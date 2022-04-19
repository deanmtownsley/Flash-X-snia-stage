!!****if* source/Simulation/SimulationForcing/incompFlow/Outlet/sim_outletData
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
!! NAME
!!
!!  sim_outletData
!!
!! SYNOPSIS
!!
!!  use sim_outletData
!!
!!***

#include "constants.h"
#include "Simulation.h"

module sim_outletData

   implicit none

   integer, save :: sim_outletFlag(LOW:HIGH, MDIM)

   real, save :: sim_outletSink
   real, save :: sim_outletBuffer
   real, save :: sim_outletGrowthRate

   real, save, dimension(LOW:HIGH, MDIM) :: sim_QOut, sim_QOutLiq, sim_QOutGas
   real, save, dimension(LOW:HIGH, MDIM) :: sim_QAux, sim_QAuxLiq, sim_QAuxGas

   real, save, dimension(LOW:HIGH, MDIM) :: sim_volOut, sim_volOutLiq, sim_volOutGas
   real, save, dimension(LOW:HIGH, MDIM) :: sim_volAux, sim_volAuxLiq, sim_volAuxGas

end module sim_outletData
