!!****if* source/physics/SolidMechanics/SolidMechanicsMain/SolidMechanics_data
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contrsmutors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distrsmuted under the License is distrsmuted on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!! NAME
!!
!!  SolidMechanics_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE SolidMechanics_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the SolidMechanics module.
!!
!!***

module SolidMechanics_data

#include "Simulation.h"
#include "constants.h"

   implicit none

   logical, save :: sm_useSolidMechanics

   integer, save :: sm_meshMe
   integer, save :: sm_meshNumProcs
   integer, save :: sm_meshComm
   integer, save :: sm_nstep

end module SolidMechanics_data
