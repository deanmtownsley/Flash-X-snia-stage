!!***if* source/Simulation/SimulationForcing/incompFlow/Outlet/phaseAvg/sim_outletVelBlkPhased
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

subroutine sim_outletVelBlk2dPhased(u, v, ru, rv, phi, xcell, ycell, &
                                    boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, &
                                    inletFlag, inletBuffer, inletGrowthRate, &
                                    outletFlag, outletBuffer, outletGrowthRate, &
                                    outletIsLiq, outletIsGas, &
                                    volAuxLiq, volAuxGas, QAuxLiq, QAuxGas, QMeanLiq, QMeanGas, &
                                    xMin, xMax, yMin, yMax, gravX, gravY)

   implicit none
   real, dimension(:, :, :), intent(in)    :: u, v
   real, dimension(:, :, :), intent(inout) :: ru, rv
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:), intent(in)          :: xcell, ycell
   real, dimension(:, :), intent(in)       :: boundBox
   real, intent(in)                        :: dt, dx, dy
   integer, intent(in)                     :: ix1, ix2, jy1, jy2
   integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag, inletFlag
   real, intent(in) :: inletBuffer, inletGrowthRate, outletBuffer, outletGrowthRate
   integer, intent(in) :: outletIsLiq, outletIsGas
   real, intent(inout) :: QAuxLiq(MDIM), QAuxGas(MDIM)
   real, intent(inout) :: volAuxLiq(MDIM), volAuxGas(MDIM)
   real, intent(in) :: QMeanLiq(MDIM), QMeanGas(MDIM)
   real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY

   integer :: i, j, k
   real :: xi, yi
   real :: uplus, umins, vplus, vmins
   real :: xforce, yforce
   real :: cellphi, cellvol
   integer :: iGas, iLiq

   cellvol = (dx/(xMax - xMin))*(dy/(yMax - yMin))

   k = 1

   do j = jy1, jy2
      do i = ix1, ix2 + 1
         xi = xcell(i)
         yi = ycell(j)

         cellphi = (phi(i, j, k) + phi(i - 1, j, k))*.5

         iLiq = (1 - int(sign(1., cellphi)))/2
         iGas = (1 + int(sign(1., cellphi)))/2

         xforce = iLiq*outletIsLiq*(QMeanLiq(IAXIS) - u(i, j, k))/dt + &
                  iGas*outletIsGas*(QMeanGas(IAXIS) - u(i, j, k))/dt
         !
         yforce = -iLiq*outletIsLiq*u(i, j, k)/dt - iGas*outletIsGas*u(i, j, k)/dt

         ru(i, j, k) = ru(i, j, k) &
                       + outletFlag(HIGH, IAXIS)*xforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                       + outletFlag(HIGH, JAXIS)*yforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer)))

         QAuxLiq(IAXIS) = QAuxLiq(IAXIS) + iLiq*outletFlag(HIGH, IAXIS)*u(i, j, k)*cellvol
         volAuxLiq(IAXIS) = volAuxLiq(IAXIS) + iLiq*outletFlag(HIGH, IAXIS)*cellvol

         QAuxGas(IAXIS) = QAuxGas(IAXIS) + iGas*outletFlag(HIGH, IAXIS)*u(i, j, k)*cellvol
         volAuxGas(IAXIS) = volAuxGas(IAXIS) + iGas*outletFlag(HIGH, IAXIS)*cellvol

      end do
   end do

   do j = jy1, jy2 + 1
      do i = ix1, ix2
         xi = xcell(i)
         yi = ycell(j)

         cellphi = (phi(i, j, k) + phi(i, j - 1, k))*.5

         iLiq = (1 - int(sign(1., cellphi)))/2
         iGas = (1 + int(sign(1., cellphi)))/2

         xforce = -iLiq*outletIsLiq*v(i, j, k)/dt - iGas*outletIsGas*v(i, j, k)/dt
         !
         yforce = iLiq*outletIsLiq*(QMeanLiq(JAXIS) - v(i, j, k))/dt + &
                  iGas*outletIsGas*(QMeanGas(JAXIS) - v(i, j, k))/dt

         rv(i, j, k) = rv(i, j, k) &
                       + outletFlag(HIGH, IAXIS)*xforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                       + outletFlag(HIGH, JAXIS)*yforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer)))

         QAuxLiq(JAXIS) = QAuxLiq(JAXIS) + iLiq*outletFlag(HIGH, JAXIS)*v(i, j, k)*cellvol
         volAuxLiq(JAXIS) = volAuxLiq(JAXIS) + iLiq*outletFlag(HIGH, JAXIS)*cellvol

         QAuxGas(JAXIS) = QAuxGas(JAXIS) + iGas*outletFlag(HIGH, JAXIS)*v(i, j, k)*cellvol
         volAuxGas(JAXIS) = volAuxGas(JAXIS) + iGas*outletFlag(HIGH, JAXIS)*cellvol

      end do
   end do

end subroutine sim_outletVelBlk2dPhased

subroutine sim_outletVelBlk3dPhased(u, v, w, ru, rv, rw, phi, xcell, ycell, zcell, &
                                    boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                    inletFlag, inletBuffer, inletGrowthRate, &
                                    outletFlag, outletBuffer, outletGrowthRate, &
                                    outletIsLiq, outletIsGas, &
                                    volAuxLiq, volAuxGas, QAuxLiq, QAuxGas, QMeanLiq, QMeanGas, &
                                    xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ)

   implicit none
   real, dimension(:, :, :), intent(in)  :: u, v, w
   real, dimension(:, :, :), intent(inout) :: ru, rv, rw
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:), intent(in)        :: xcell, ycell, zcell
   real, dimension(:, :), intent(in)     :: boundBox
   real, intent(in)                      :: dt, dx, dy, dz
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
   integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag, inletFlag
   real, intent(in) :: inletBuffer, inletGrowthRate, outletBuffer, outletGrowthRate
   integer, intent(in) :: outletIsLiq, outletIsGas
   real, intent(inout) :: QAuxLiq(MDIM), QAuxGas(MDIM)
   real, intent(inout) :: volAuxLiq(MDIM), volAuxGas(MDIM)
   real, intent(in) :: QMeanLiq(MDIM), QMeanGas(MDIM)
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ

   integer :: i, j, k
   real    :: xi, yi, zi
   real    :: uplus, umins, vplus, vmins, wplus, wmins
   real    :: xforce, yforce, zforce
   real    :: cellphi, cellvol
   integer :: iLiq, iGas

   cellvol = (dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2 + 1
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            cellphi = (phi(i, j, k) + phi(i - 1, j, k))*.5

            iLiq = (1 - int(sign(1., cellphi)))/2
            iGas = (1 + int(sign(1., cellphi)))/2

            xforce = iLiq*outletIsLiq*(QMeanLiq(IAXIS) - u(i, j, k))/dt + &
                     iGas*outletIsGas*(QMeanGas(IAXIS) - u(i, j, k))/dt
            !
            yforce = -iLiq*outletIsLiq*u(i, j, k)/dt - iGas*outletIsGas*u(i, j, k)/dt
            zforce = -iLiq*outletIsLiq*u(i, j, k)/dt - iGas*outletIsGas*u(i, j, k)/dt

            ru(i, j, k) = ru(i, j, k) &
                          + outletFlag(HIGH, IAXIS)*xforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                          + outletFlag(HIGH, JAXIS)*yforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer))) &
                          + outletFlag(HIGH, KAXIS)*zforce*(2/(1 + exp(-outletGrowthRate*(zi - zMax)/outletBuffer)))

            QAuxLiq(IAXIS) = QAuxLiq(IAXIS) + iLiq*outletFlag(HIGH, IAXIS)*u(i, j, k)*cellvol
            volAuxLiq(IAXIS) = volAuxLiq(IAXIS) + iLiq*outletFlag(HIGH, IAXIS)*cellvol

            QAuxGas(IAXIS) = QAuxGas(IAXIS) + iGas*outletFlag(HIGH, IAXIS)*u(i, j, k)*cellvol
            volAuxGas(IAXIS) = volAuxGas(IAXIS) + iGas*outletFlag(HIGH, IAXIS)*cellvol

         end do
      end do
   end do

   do k = kz1, kz2
      do j = jy1, jy2 + 1
         do i = ix1, ix2
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            cellphi = (phi(i, j, k) + phi(i, j - 1, k))*.5

            iLiq = (1 - int(sign(1., cellphi)))/2
            iGas = (1 + int(sign(1., cellphi)))/2

            xforce = -iLiq*outletIsLiq*v(i, j, k)/dt - iGas*outletIsGas*v(i, j, k)/dt
            !
            yforce = iLiq*outletIsLiq*(QMeanLiq(JAXIS) - v(i, j, k))/dt + &
                     iGas*outletIsGas*(QMeanGas(JAXIS) - v(i, j, k))/dt
            !
            zforce = -iLiq*outletIsLiq*v(i, j, k)/dt - iGas*outletIsGas*v(i, j, k)/dt

            rv(i, j, k) = rv(i, j, k) &
                          + outletFlag(HIGH, IAXIS)*xforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                          + outletFlag(HIGH, JAXIS)*yforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer))) &
                          + outletFlag(HIGH, KAXIS)*zforce*(2/(1 + exp(-outletGrowthRate*(zi - zMax)/outletBuffer)))

            QAuxLiq(JAXIS) = QAuxLiq(JAXIS) + iLiq*outletFlag(HIGH, JAXIS)*v(i, j, k)*cellvol
            volAuxLiq(JAXIS) = volAuxLiq(JAXIS) + iLiq*outletFlag(HIGH, JAXIS)*cellvol

            QAuxGas(JAXIS) = QAuxGas(JAXIS) + iGas*outletFlag(HIGH, JAXIS)*v(i, j, k)*cellvol
            volAuxGas(JAXIS) = volAuxGas(JAXIS) + iGas*outletFlag(HIGH, JAXIS)*cellvol

         end do
      end do
   end do

   do k = kz1, kz2 + 1
      do j = jy1, jy2
         do i = ix1, ix2
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            cellphi = (phi(i, j, k) + phi(i, j, k - 1))*.5

            iLiq = (1 - int(sign(1., cellphi)))/2
            iGas = (1 + int(sign(1., cellphi)))/2

            xforce = -iLiq*outletIsLiq*w(i, j, k)/dt - iGas*outletIsGas*w(i, j, k)/dt
            yforce = -iLiq*outletIsLiq*w(i, j, k)/dt - iGas*outletIsGas*w(i, j, k)/dt
            !
            zforce = iLiq*outletIsLiq*(QMeanLiq(KAXIS) - w(i, j, k))/dt + &
                     iGas*outletIsGas*(QMeanGas(KAXIS) - w(i, j, k))/dt

            rw(i, j, k) = rw(i, j, k) &
                          + outletFlag(HIGH, IAXIS)*xforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                          + outletFlag(HIGH, JAXIS)*yforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer))) &
                          + outletFlag(HIGH, KAXIS)*zforce*(2/(1 + exp(-outletGrowthRate*(zi - zMax)/outletBuffer)))

            QAuxLiq(KAXIS) = QAuxLiq(KAXIS) + iLiq*outletFlag(HIGH, KAXIS)*w(i, j, k)*cellvol
            volAuxLiq(KAXIS) = volAuxLiq(KAXIS) + iLiq*outletFlag(HIGH, KAXIS)*cellvol

            QAuxGas(KAXIS) = QAuxGas(KAXIS) + iGas*outletFlag(HIGH, KAXIS)*w(i, j, k)*cellvol
            volAuxGas(KAXIS) = volAuxGas(KAXIS) + iGas*outletFlag(HIGH, KAXIS)*cellvol

         end do
      end do
   end do

end subroutine sim_outletVelBlk3dPhased
