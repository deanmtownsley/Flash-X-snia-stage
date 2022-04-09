!!***if* source/Simulation/SimulationForcing/incompFlow/Outlet/sim_outletVelBlk
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
                              boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, &
                              inletFlag, inletBuffer, inletGrowthRate, &
                              outletFlag, outletBuffer, outletGrowthRate, &
                              volAux, QAux, QMean, &
                              xMin, xMax, yMin, yMax, gravX, gravY)

   implicit none
   real, dimension(:, :, :), intent(in)    :: u, v
   real, dimension(:, :, :), intent(inout) :: ru, rv
   real, dimension(:), intent(in)          :: xcell, ycell
   real, dimension(:, :), intent(in)       :: boundBox
   real, intent(in)                        :: dt, dx, dy
   integer, intent(in)                     :: ix1, ix2, jy1, jy2
   integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag, inletFlag
   real, intent(in) :: inletBuffer, inletGrowthRate, outletBuffer, outletGrowthRate
   real, intent(inout) :: QAux(MDIM), volAux(MDIM)
   real, intent(in) :: QMean(MDIM)
   real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY

   integer :: i, j, k
   real    :: xi, yi
   real    :: uplus, umins, vplus, vmins
   real    :: xforce, yforce, cellvol

   cellvol = (dx/(xMax - xMin))*(dy/(yMax - yMin))

   k = 1

   do j = jy1, jy2
      do i = ix1, ix2 + 1
         xi = xcell(i)
         yi = ycell(j)

         xforce = (QMean(IAXIS) - u(i, j, k))/dt
         yforce = -u(i, j, k)/dt

         ru(i, j, k) = ru(i, j, k) &
                       + outletFlag(HIGH, IAXIS)*xforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                       + outletFlag(HIGH, JAXIS)*yforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer)))

         QAux(IAXIS) = QAux(IAXIS) + &
                       outletFlag(HIGH, IAXIS)*u(i, j, k)*cellvol

         volAux(IAXIS) = volAux(IAXIS) + &
                         outletFlag(HIGH, IAXIS)*cellvol
      end do
   end do

   do j = jy1, jy2 + 1
      do i = ix1, ix2
         xi = xcell(i)
         yi = ycell(j)

         xforce = -v(i, j, k)/dt
         yforce = (QMean(JAXIS) - v(i, j, k))/dt

         rv(i, j, k) = rv(i, j, k) &
                       + outletFlag(HIGH, IAXIS)*xforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                       + outletFlag(HIGH, JAXIS)*yforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer)))

         QAux(JAXIS) = QAux(JAXIS) + &
                       outletFlag(HIGH, JAXIS)*v(i, j, k)*cellvol

         volAux(JAXIS) = volAux(JAXIS) + &
                         outletFlag(HIGH, JAXIS)*cellvol

      end do
   end do

end subroutine sim_outletVelBlk2d

subroutine sim_outletVelBlk3d(u, v, w, ru, rv, rw, xcell, ycell, zcell, &
                              boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                              inletFlag, inletBuffer, inletGrowthRate, &
                              outletFlag, outletBuffer, outletGrowthRate, &
                              volAux, QAux, QMean, &
                              xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ)

   implicit none
   real, dimension(:, :, :), intent(in)  :: u, v, w
   real, dimension(:, :, :), intent(inout) :: ru, rv, rw
   real, dimension(:), intent(in)        :: xcell, ycell, zcell
   real, dimension(:, :), intent(in)     :: boundBox
   real, intent(in)                      :: dt, dx, dy, dz
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
   integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag, inletFlag
   real, intent(in) :: inletBuffer, inletGrowthRate, outletBuffer, outletGrowthRate
   real, intent(inout) :: QAux(MDIM), volAux(MDIM)
   real, intent(in) :: QMean(MDIM)
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ

   integer :: i, j, k
   real    :: xi, yi, zi
   real    :: uplus, umins, vplus, vmins, wplus, wmins
   real    :: xforce, yforce, zforce, cellvol

   cellvol = (dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2 + 1
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            xforce = (QMean(IAXIS) - u(i, j, k))/dt
            yforce = -u(i, j, k)/dt
            zforce = -u(i, j, k)/dt

            ru(i, j, k) = ru(i, j, k) &
                          + outletFlag(HIGH, IAXIS)*xforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                          + outletFlag(HIGH, JAXIS)*yforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer))) &
                          + outletFlag(HIGH, KAXIS)*zforce*(2/(1 + exp(-outletGrowthRate*(zi - zMax)/outletBuffer)))

            QAux(IAXIS) = QAux(IAXIS) + &
                          outletFlag(HIGH, IAXIS)*u(i, j, k)*cellvol

            volAux(IAXIS) = volAux(IAXIS) + &
                            outletFlag(HIGH, IAXIS)*cellvol

         end do
      end do
   end do

   do k = kz1, kz2
      do j = jy1, jy2 + 1
         do i = ix1, ix2
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            xforce = -v(i, j, k)/dt
            yforce = (QMean(JAXIS) - v(i, j, k))/dt
            zforce = -v(i, j, k)/dt

            rv(i, j, k) = rv(i, j, k) &
                          + outletFlag(HIGH, IAXIS)*xforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                          + outletFlag(HIGH, JAXIS)*yforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer))) &
                          + outletFlag(HIGH, KAXIS)*zforce*(2/(1 + exp(-outletGrowthRate*(zi - zMax)/outletBuffer)))

            QAux(JAXIS) = QAux(JAXIS) + &
                          outletFlag(HIGH, JAXIS)*v(i, j, k)*cellvol

            volAux(JAXIS) = volAux(JAXIS) + &
                            outletFlag(HIGH, JAXIS)*cellvol

         end do
      end do
   end do

   do k = kz1, kz2 + 1
      do j = jy1, jy2
         do i = ix1, ix2
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            xforce = -w(i, j, k)/dt
            yforce = -w(i, j, k)/dt
            zforce = (QMean(KAXIS) - w(i, j, k))/dt

            rw(i, j, k) = rw(i, j, k) &
                          + outletFlag(HIGH, IAXIS)*xforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                          + outletFlag(HIGH, JAXIS)*yforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer))) &
                          + outletFlag(HIGH, KAXIS)*zforce*(2/(1 + exp(-outletGrowthRate*(zi - zMax)/outletBuffer)))

            QAux(KAXIS) = QAux(KAXIS) + &
                          outletFlag(HIGH, KAXIS)*w(i, j, k)*cellvol

            volAux(KAXIS) = volAux(KAXIS) + &
                            outletFlag(HIGH, KAXIS)*cellvol

         end do
      end do
   end do

end subroutine sim_outletVelBlk3d
