!!***if* source/Simulation/SimulationForcing/incompFlow/Outflow/sim_outflowLSDampingBlk
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

subroutine sim_outflowLSDampingBlk2d(pfrc, phi, xcell, ycell, boundBox, &
                                     dt, dx, dy, ix1, ix2, jy1, jy2, &
                                     domainBC, outflowSink, outflowBuffer, outflowGrowthRate, &
                                     xMin, xMax, yMin, yMax)

   implicit none

   real, dimension(:, :, :), intent(inout) :: pfrc
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:), intent(in) :: xcell, ycell
   real, dimension(:, :), intent(in) :: boundBox
   real, intent(in) :: dt, dx, dy
   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer, dimension(2, MDIM), intent(in) :: domainBC
   real, intent(in) :: outflowSink, outflowBuffer, outflowGrowthRate
   real, intent(in) :: xMin, xMax, yMin, yMax

   integer :: i, j, k
   real    :: xi, yi

   k = 1

   do j = jy1, jy2
      do i = ix1, ix2
         xi = xcell(i)
         yi = ycell(j)

         if (domainBC(HIGH, JAXIS) == NEUMANN_INS .or. &
             domainBC(HIGH, JAXIS) == OUTFLOW_INS .or. &
             domainBC(HIGH, JAXIS) == EXTRAP_INS) then
            pfrc(i, j, k) = pfrc(i, j, k) - outflowSink* &
                            (2/(1 + exp(-outflowGrowthRate*(yi - yMax)/outflowBuffer)))

         else if (domainBC(HIGH, IAXIS) == NEUMANN_INS .or. &
                  domainBC(HIGH, IAXIS) == OUTFLOW_INS .or. &
                  domainBC(HIGH, IAXIS) == EXTRAP_INS) then
            pfrc(i, j, k) = pfrc(i, j, k) - outflowSink* &
                            (2/(1 + exp(-outflowGrowthRate*(xi - yMax)/outflowBuffer)))

         end if

      end do
   end do

end subroutine sim_outflowLSDampingBlk2d

subroutine sim_outflowLSDampingBlk3d(pfrc, phi, xcell, ycell, zcell, boundBox, &
                                     dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                     domainBC, outflowSink, outflowBuffer, outflowGrowthRate, &
                                     xMin, xMax, yMin, yMax, zMin, zMax)

   implicit none

   real, dimension(:, :, :), intent(inout) :: pfrc
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:), intent(in) :: xcell, ycell, zcell
   real, dimension(:, :), intent(in) :: boundBox
   real, intent(in) :: dt, dx, dy, dz
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
   integer, dimension(2, MDIM), intent(in) :: domainBC
   real, intent(in) :: outflowSink, outflowBuffer, outflowGrowthRate
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax

   integer :: i, j, k
   real    :: xi, yi, zi

   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            if (domainBC(HIGH, JAXIS) == NEUMANN_INS .or. &
                domainBC(HIGH, JAXIS) == OUTFLOW_INS .or. &
                domainBC(HIGH, JAXIS) == EXTRAP_INS) then
               pfrc(i, j, k) = pfrc(i, j, k) - outflowSink* &
                               (2/(1 + exp(-outflowGrowthRate*(yi - yMax)/outflowBuffer)))

            else if (domainBC(HIGH, IAXIS) == NEUMANN_INS .or. &
                     domainBC(HIGH, IAXIS) == OUTFLOW_INS .or. &
                     domainBC(HIGH, IAXIS) == EXTRAP_INS) then
               pfrc(i, j, k) = pfrc(i, j, k) - outflowSink* &
                               (2/(1 + exp(-outflowGrowthRate*(xi - xMax)/outflowBuffer)))

            else if (domainBC(HIGH, KAXIS) == NEUMANN_INS .or. &
                     domainBC(HIGH, KAXIS) == OUTFLOW_INS .or. &
                     domainBC(HIGH, KAXIS) == EXTRAP_INS) then
               pfrc(i, j, k) = pfrc(i, j, k) - outflowSink* &
                               (2/(1 + exp(-outflowGrowthRate*(zi - zMax)/outflowBuffer)))

            end if

         end do
      end do
   end do

end subroutine sim_outflowLSDampingBlk3d
