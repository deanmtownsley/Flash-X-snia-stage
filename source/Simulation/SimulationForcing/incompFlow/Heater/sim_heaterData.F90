!!****if* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterData
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
!!  sim_heaterData
!!
!! SYNOPSIS
!!
!!  use sim_heaterData
!!
!!***

#include "constants.h"
#include "Simulation.h"

module sim_heaterData

   implicit none

   real, save        :: sim_nucSeedRadius
   integer, save     :: sim_numHeaters
   character(len=20), save :: sim_heaterName
   logical, save :: sim_heaterShowInfo

   type sim_heaterType

      real    :: xMin
      real    :: xMax
      real    :: zMin
      real    :: zMax
      real    :: yMin
      real    :: yMax

      real    :: wallTemp
      real    :: nucWaitTime
      real    :: advAngle
      real    :: rcdAngle
      real    :: seedRadius
      real    :: seedHeight
      real    :: velContact

      integer :: numSites

      real, dimension(:), allocatable :: xSite
      real, dimension(:), allocatable :: zSite
      real, dimension(:), allocatable :: ySite
      real, dimension(:), allocatable :: siteRadii
      real, dimension(:), allocatable :: siteTimeStamp
      logical, dimension(:), allocatable :: siteIsAttachedCurr
      logical, dimension(:), allocatable :: siteIsAttachedPrev

   end type sim_heaterType

   type(sim_heaterType), save, dimension(:), pointer :: sim_heaterInfo

end module sim_heaterData
