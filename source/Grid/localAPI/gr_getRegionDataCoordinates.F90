!!****if* source/Grid/localAPI/gr_getRegionDataCoordinates.F90
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
!!  gr_getRegionDataCoordinates
!!
!! SYNOPSIS
!!  call Grid_bcApplyToRegionSpecialized(integer(IN)  :: level, 
!!                                       integer(IN)  :: gridDataStruct,
!!                                       integer(IN)  :: axis,
!!                                       integer(IN)  :: axis2,
!!                                       integer(IN)  :: axis3,
!!                                       integer(IN)  :: regionSize(REGION_DIM),
!!                                       integer(IN)  :: endPoints(LOW:HIGH, MDIM),
!!                                       real(OUT)    :: coordinates(regionSize(BC_DIR),
!!                                                                   regionSize(SECOND_DIR), 
!!                                                                   regionSize(THIRD_DIR),
!!                                                                   1:MDIM)) 
!!
!! DESCRIPTION
!!  This is a helper routine for those versions of Grid_bcApplyToRegionSpecialized
!!  that require knowledge of the physical coordinates of all mesh elements in
!!  the given regionData on which they are to operate.
!!
!! ARGUMENTS
!!  See ARGUMENTS section of documentation for Grid_bcApplyToRegion.
!!
!!***

#include "constants.h"

subroutine gr_getRegionDataCoordinates(level, gridDataStruct, &
                                       axis, axis2, axis3, &
                                       regionSize, endPoints, &
                                       coordinates)
  implicit none

  integer, intent(IN)  :: level
  integer, intent(IN)  :: gridDataStruct
  integer, intent(IN)  :: axis
  integer, intent(IN)  :: axis2
  integer, intent(IN)  :: axis3
  integer, intent(IN)  :: regionSize(REGION_DIM)
  integer, intent(IN)  :: endPoints(LOW:HIGH, 1:MDIM)
  real,    intent(OUT) :: coordinates(regionSize(BC_DIR),     &
                                      regionSize(SECOND_DIR), &
                                      regionSize(THIRD_DIR),  &
                                      1:MDIM)
  
  coordinates(:, :, :, :) = 0.0
end subroutine gr_getRegionDataCoordinates

