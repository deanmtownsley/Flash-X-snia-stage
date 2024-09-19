!!****if* source/physics/sourceTerms/Heater/HeaterMain/Heater_type
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
!!  Heater_type
!!
!! SYNOPSIS
!!
!!  use Heater_type
!!
!!***

#include "constants.h"
#include "Simulation.h"

module Heater_type

   use iso_c_binding
   implicit none

   type Heater_type_t

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

      real    :: C3 = 0.
      real    :: C2 = 0.
      real    :: C1 = 0.
      real    :: C0 = 0.
      real    :: tbl_thickness = 0.2
      integer :: non_uniform_temp_flag = 0
      integer :: heat_flux_flag = 0
      real    :: nd_heat_flux = 0.

      integer :: numSitesAll, numSitesProc
      integer :: numSitesBlk(MAXBLOCKS), siteMapOnProc(MAXBLOCKS, HTR_MAX_NUMSITES)

      real, dimension(:), allocatable :: xSiteInit, ySiteInit, zSiteInit, radiusInit
      real, dimension(HTR_MAX_NUMSITES) :: xSiteProc, zSiteProc, ySiteProc, siteTimeStamp
      logical, dimension(HTR_MAX_NUMSITES) :: siteIsAttachedCurr, siteIsAttachedPrev

#ifdef HEATER_ANN_SEARCH
      type(c_ptr) :: kdTree = c_null_ptr
#endif
      integer :: dims

   end type Heater_type_t

end module Heater_type
