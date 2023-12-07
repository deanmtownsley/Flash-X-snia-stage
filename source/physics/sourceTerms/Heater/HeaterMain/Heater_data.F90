!!****if* source/physics/sourceTerms/Heater/HeaterMain/Heater_data
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
!!  Heater_data
!!
!! SYNOPSIS
!!
!!  use Heater_data
!!
!!***

#include "constants.h"
#include "Simulation.h"

module Heater_data

   use Heater_type, ONLY: Heater_type_t

   implicit none

   integer, save     :: htr_meshMe
   real, save        :: htr_nucSeedRadius
   integer, save     :: htr_numHeaters
   character(len=20), save :: htr_heaterName
   logical, save :: htr_showInfo

   type(Heater_type_t), save, dimension(:), pointer :: htr_heaterInfo

#ifdef HEATER_ANN_SEARCH
   integer, save, dimension(:), allocatable :: htr_annIdx
   integer, save :: htr_annQueries
#endif

end module Heater_data
