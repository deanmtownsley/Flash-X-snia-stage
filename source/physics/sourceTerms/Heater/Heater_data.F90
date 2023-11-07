!!****if* source/Simulation/SimulationForcing/incompFlow/Heater/Heater_Data
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
!!  Heater_Data
!!
!! SYNOPSIS
!!
!!  use Heater_Data
!!
!!***

#include "constants.h"
#include "Simulation.h"

module Heater_Data

   use iso_c_binding
   implicit none

   real, save        :: sim_nucSeedRadius
   integer, save     :: sim_numHeaters
   character(len=20), save :: Heater_Name
   logical, save :: Heater_ShowInfo

   type Heater_Type

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

      integer :: numSitesAll, numSitesProc
      integer :: numSitesBlk(MAXBLOCKS), siteMapOnProc(MAXBLOCKS, SIM_MAX_NUMSITES)

      real, dimension(:), allocatable :: xSiteInit, ySiteInit, zSiteInit, radiusInit
      real, dimension(SIM_MAX_NUMSITES) :: xSiteProc, zSiteProc, ySiteProc, siteTimeStamp
      logical, dimension(SIM_MAX_NUMSITES) :: siteIsAttachedCurr, siteIsAttachedPrev

#ifdef SIM_HEATER_ANN_SEARCH
      type(c_ptr) :: kdTree = c_null_ptr
#endif
      integer :: dims

   end type Heater_Type

   type(Heater_Type), save, dimension(:), pointer :: Heater_Info

#ifdef SIM_HEATER_ANN_SEARCH
   integer, save, dimension(:), allocatable :: Heater_AnnIdx
   integer, save :: Heater_AnnQueries
#endif

end module Heater_Data
