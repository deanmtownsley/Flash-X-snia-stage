!!***if* source/Simulation/SimulationForcing/incompFlow/Outlet/sim_outletLSDampingBlk
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

subroutine sim_outletLSDampingBlk2d(pfrc, phi, xcell, ycell, boundBox, &
                                     dt, dx, dy, ix1, ix2, jy1, jy2, &
                                     domainBC, outletSink, outletBuffer, outletGrowthRate, &
                                     xMin, xMax, yMin, yMax)

   implicit none

   real, dimension(:, :, :), intent(inout) :: pfrc
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:), intent(in) :: xcell, ycell
   real, dimension(:, :), intent(in) :: boundBox
   real, intent(in) :: dt, dx, dy
   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer, dimension(2, MDIM), intent(in) :: domainBC
   real, intent(in) :: outletSink, outletBuffer, outletGrowthRate
   real, intent(in) :: xMin, xMax, yMin, yMax

   integer :: i, j, k
   real    :: xi, yi

end subroutine sim_outletLSDampingBlk2d

subroutine sim_outletLSDampingBlk3d(pfrc, phi, xcell, ycell, zcell, boundBox, &
                                     dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                     domainBC, outletSink, outletBuffer, outletGrowthRate, &
                                     xMin, xMax, yMin, yMax, zMin, zMax)

   implicit none

   real, dimension(:, :, :), intent(inout) :: pfrc
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:), intent(in) :: xcell, ycell, zcell
   real, dimension(:, :), intent(in) :: boundBox
   real, intent(in) :: dt, dx, dy, dz
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
   integer, dimension(2, MDIM), intent(in) :: domainBC
   real, intent(in) :: outletSink, outletBuffer, outletGrowthRate
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax

   integer :: i, j, k
   real    :: xi, yi, zi

end subroutine sim_outletLSDampingBlk3d
