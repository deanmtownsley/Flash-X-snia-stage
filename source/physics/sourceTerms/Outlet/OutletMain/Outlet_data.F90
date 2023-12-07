!!****if* source/physics/sourceTerms/Outlet/OutletMain/Outlet_data
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
!! NAME
!!
!!  Outlet_data
!!
!! SYNOPSIS
!!
!!  use Outlet_data
!!
!!***

#include "constants.h"
#include "Simulation.h"

module Outlet_data

   implicit none

   integer, save :: out_meshMe
   real, save :: out_xMin, out_xMax, out_yMin, out_yMax, out_zMin, out_zMax

   integer, save :: out_flag(LOW:HIGH, MDIM)

   real, save :: out_sink
   real, save :: out_buffer
   real, save :: out_growthRate
   real, save :: out_velRefScale

   real, save, dimension(LOW:HIGH, MDIM) :: out_QOut, out_QOutLiq, out_QOutGas
   real, save, dimension(LOW:HIGH, MDIM) :: out_QAux, out_QAuxLiq, out_QAuxGas

   real, save, dimension(LOW:HIGH, MDIM) :: out_volOut, out_volOutLiq, out_volOutGas
   real, save, dimension(LOW:HIGH, MDIM) :: out_volAux, out_volAuxLiq, out_volAuxGas

end module Outlet_data
