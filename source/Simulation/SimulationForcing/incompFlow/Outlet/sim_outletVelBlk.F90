!!***if* source/Simulation/SimulationForcing/incompFlow/Outlet/velOutBlk
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

subroutine sim_outletVelBlk2d(u, v, ru, rv, xcell, ycell, &
                              boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, inletFlag, &
                              outletFlag, velOutAux, velOut, outletBuffer, outletGrowthRate, &
                              xMin, xMax, yMin, yMax, gravX, gravY)

   implicit none
   real, dimension(:, :, :), intent(in)    :: u, v
   real, dimension(:, :, :), intent(inout) :: ru, rv
   real, dimension(:), intent(in)          :: xcell, ycell
   real, dimension(:, :), intent(in)       :: boundBox
   real, intent(in)                        :: dt, dx, dy
   integer, intent(in)                     :: ix1, ix2, jy1, jy2
   integer, dimension(2, MDIM), intent(in) :: outletFlag, inletFlag
   real, intent(inout) :: velOutAux(2, MDIM)
   real, intent(in) :: velOut(2, MDIM), outletBuffer, outletGrowthRate
   real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY

   integer :: i, j, k
   real    :: xi, yi
   real    :: uplus, umins, vplus, vmins
   real    :: xforce, yforce

   k = 1

   do j = jy1, jy2
      do i = ix1, ix2 + 1
         xi = xcell(i)
         yi = ycell(j)

         xforce = (velOut(HIGH, IAXIS) - u(i, j, k))/dt
         yforce = (velOut(HIGH, IAXIS) - u(i, j, k))/dt

         ru(i, j, k) = ru(i, j, k) &
                       + outletFlag(HIGH, IAXIS)*xforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                       + outletFlag(HIGH, JAXIS)*yforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer)))

         velOutAux(HIGH, IAXIS) = velOutAux(HIGH, IAXIS) + &
                                  outletFlag(HIGH, IAXIS)*u(i, j, k)*(dx/(xMax - xMin))*(dy/(yMax - yMin))

      end do
   end do

   do j = jy1, jy2 + 1
      do i = ix1, ix2
         xi = xcell(i)
         yi = ycell(j)

         xforce = (velOut(HIGH, JAXIS) - v(i, j, k))/dt
         yforce = (velOut(HIGH, JAXIS) - v(i, j, k))/dt

         rv(i, j, k) = rv(i, j, k) &
                       + outletFlag(HIGH, IAXIS)*xforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                       + outletFlag(HIGH, JAXIS)*yforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer)))

         velOutAux(HIGH, JAXIS) = velOutAux(HIGH, JAXIS) + &
                                  outletFlag(HIGH, JAXIS)*v(i, j, k)*(dx/(xMax - xMin))*(dy/(yMax - yMin))

      end do
   end do

end subroutine sim_outletVelBlk2d

subroutine sim_outletVelBlk3d(u, v, w, ru, rv, rw, xcell, ycell, zcell, &
                              boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, inletFlag, &
                              outletFlag, velOutAux, velOut, outletBuffer, outletGrowthRate, &
                              xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ)

   implicit none
   real, dimension(:, :, :), intent(in)  :: u, v, w
   real, dimension(:, :, :), intent(inout) :: ru, rv, rw
   real, dimension(:), intent(in)        :: xcell, ycell, zcell
   real, dimension(:, :), intent(in)     :: boundBox
   real, intent(in)                      :: dt, dx, dy, dz
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
   integer, dimension(2, MDIM), intent(in) :: outletFlag, inletFlag
   real, intent(inout) :: velOutAux(2, MDIM)
   real, intent(in) :: velOut(2, MDIM), outletBuffer, outletGrowthRate
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ

   integer :: i, j, k
   real    :: xi, yi, zi
   real    :: uplus, umins, vplus, vmins, wplus, wmins

   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2 + 1
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

         end do
      end do
   end do

   do k = kz1, kz2
      do j = jy1, jy2 + 1
         do i = ix1, ix2
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

         end do
      end do
   end do

   do k = kz1, kz2 + 1
      do j = jy1, jy2
         do i = ix1, ix2
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

         end do
      end do
   end do

end subroutine sim_outletVelBlk3d
