!!***if* source/Simulation/SimulationForcing/incompFlow/Outlet/sim_outletLSDamping
!!
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
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine sim_outletLSDamping(pfrc, phi, xcenter, ycenter, zcenter, boundBox, &
                                 dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                 outletFlag, outletSink, outletBuffer, outletGrowthRate, &
                                 xMin, xMax, yMin, yMax, zMin, zMax)

   implicit none

   real, dimension(:, :, :), intent(inout) :: pfrc
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:), intent(in) :: xcenter, ycenter, zcenter
   real, dimension(:, :), intent(in) :: boundBox
   real, intent(in) :: dt, dx, dy, dz
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
   integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag
   real, intent(in) :: outletSink, outletBuffer, outletGrowthRate
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax

   integer :: i, j, k
   real    :: xi, yi, zi

end subroutine sim_outletLSDamping
