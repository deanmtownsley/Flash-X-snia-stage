!!***if* source/Simulation/SimulationForcing/incompFlow/Outlet/phaseForcing/sim_outletVelBlkPhased
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

subroutine sim_outletVelBlk2dPhased(u, v, ru, rv, phi, xcenter, ycenter, &
                                    boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, &
                                    inletFlag, inletBuffer, inletGrowthRate, &
                                    outletFlag, outletBuffer, outletGrowthRate, &
                                    volAuxLiq, volAuxGas, QAuxLiq, QAuxGas, QOutLiq, QOutGas, &
                                    xMin, xMax, yMin, yMax, gravX, gravY)

   implicit none
   real, dimension(:, :, :), intent(in)    :: u, v
   real, dimension(:, :, :), intent(inout) :: ru, rv
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:), intent(in)          :: xcenter, ycenter
   real, dimension(:, :), intent(in)       :: boundBox
   real, intent(in)                        :: dt, dx, dy
   integer, intent(in)                     :: ix1, ix2, jy1, jy2
   integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag, inletFlag
   real, intent(in) :: inletBuffer, inletGrowthRate, outletBuffer, outletGrowthRate
   real, intent(inout) :: QAuxLiq(LOW:HIGH, MDIM), QAuxGas(LOW:HIGH, MDIM)
   real, intent(inout) :: volAuxLiq(LOW:HIGH, MDIM), volAuxGas(LOW:HIGH, MDIM)
   real, intent(in) :: QOutLiq(LOW:HIGH, MDIM), QOutGas(LOW:HIGH, MDIM)
   real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY

   integer :: i, j, k, idimn, ibound, iliq, igas, iforce(LOW:HIGH)
   real :: xcell, ycell, volcell, phicell
   real :: velforce(LOW:HIGH, MDIM), outprofile(LOW:HIGH, MDIM), velout(LOW:HIGH)
   real, parameter :: velref = 1.

   volcell = dx*dy

   k = 1

   do j = jy1, jy2
      do i = ix1, ix2 + 1
         xcell = xcenter(i) - dx/2
         ycell = ycenter(j)

         phicell = (phi(i, j, k) + phi(i - 1, j, k))*.5
         iliq = (1 - int(sign(1., phicell)))/2
         igas = (1 + int(sign(1., phicell)))/2

         outprofile(LOW, IAXIS:JAXIS) = (/2/(1 + exp(outletGrowthRate*(xcell - xMin)/outletBuffer)), &
                                          2/(1 + exp(outletGrowthRate*(ycell - yMin)/outletBuffer))/)
         !
         outprofile(HIGH, IAXIS:JAXIS) = (/2/(1 + exp(-outletGrowthRate*(xcell - xMax)/outletBuffer)), &
                                           2/(1 + exp(-outletGrowthRate*(ycell - yMax)/outletBuffer))/)

         velout(LOW) = iliq*abs(QOutLiq(LOW, IAXIS)) + igas*abs(QOutGas(LOW, IAXIS))
         velout(HIGH) = iliq*abs(QOutLiq(HIGH, IAXIS)) + igas*abs(QOutGas(HIGH, IAXIS))

         iforce(LOW) = abs(u(i, j, k)) > max(velref, velout(LOW))
         iforce(HIGH) = abs(u(i, j, k)) > max(velref, velout(HIGH))

         velforce(LOW:HIGH, IAXIS) = (/iforce(LOW)*(velout(LOW)*u(i, j, k)/(abs(u(i, j, k)) + 1e-13) - u(i, j, k))/dt, &
                                       iforce(HIGH)*(velout(HIGH)*u(i, j, k)/(abs(u(i, j, k)) + 1e-13) - u(i, j, k))/dt/)
         !
         velforce(LOW:HIGH, JAXIS) = -u(i, j, k)/dt

         do idimn = 1, NDIM
            do ibound = LOW, HIGH
               ru(i, j, k) = ru(i, j, k) + outletFlag(ibound, idimn)*velforce(ibound, idimn)*outprofile(ibound, idimn)
            end do
         end do

         do ibound = LOW, HIGH
            QAuxLiq(ibound, IAXIS) = QAuxLiq(ibound, IAXIS) + &
                                     iLiq*outletFlag(ibound, IAXIS)*u(i, j, k)*volcell*outprofile(ibound, IAXIS)

            volAuxLiq(ibound, IAXIS) = volAuxLiq(ibound, IAXIS) + &
                                       iLiq*outletFlag(ibound, IAXIS)*volcell*outprofile(ibound, IAXIS)

            QAuxGas(ibound, IAXIS) = QAuxGas(ibound, IAXIS) + &
                                     iGas*outletFlag(ibound, IAXIS)*u(i, j, k)*volcell*outprofile(ibound, IAXIS)

            volAuxGas(ibound, IAXIS) = volAuxGas(ibound, IAXIS) + &
                                       iGas*outletFlag(ibound, IAXIS)*volcell*outprofile(ibound, IAXIS)
         end do

      end do
   end do

   do j = jy1, jy2 + 1
      do i = ix1, ix2
         xcell = xcenter(i)
         ycell = ycenter(j) - dy/2

         phicell = (phi(i, j, k) + phi(i, j - 1, k))*.5
         iliq = (1 - int(sign(1., phicell)))/2
         igas = (1 + int(sign(1., phicell)))/2

         outprofile(LOW, IAXIS:JAXIS) = (/2/(1 + exp(outletGrowthRate*(xcell - xMin)/outletBuffer)), &
                                          2/(1 + exp(outletGrowthRate*(ycell - yMin)/outletBuffer))/)
         !
         outprofile(HIGH, IAXIS:JAXIS) = (/2/(1 + exp(-outletGrowthRate*(xcell - xMax)/outletBuffer)), &
                                           2/(1 + exp(-outletGrowthRate*(ycell - yMax)/outletBuffer))/)

         velout(LOW) = iliq*abs(QOutLiq(LOW, JAXIS)) + igas*abs(QOutGas(LOW, JAXIS))
         velout(HIGH) = iliq*abs(QOutLiq(HIGH, JAXIS)) + igas*abs(QOutGas(HIGH, JAXIS))

         iforce(LOW) = abs(v(i, j, k)) > max(velref, velout(LOW))
         iforce(HIGH) = abs(v(i, j, k)) > max(velref, velout(HIGH))

         velforce(LOW:HIGH, IAXIS) = -v(i, j, k)/dt
         !
         velforce(LOW:HIGH, JAXIS) = (/iforce(LOW)*(velout(LOW)*v(i, j, k)/(abs(v(i, j, k)) + 1e-13) - v(i, j, k))/dt, &
                                       iforce(HIGH)*(velout(HIGH)*v(i, j, k)/(abs(v(i, j, k)) + 1e-13) - v(i, j, k))/dt/)

         do idimn = 1, NDIM
            do ibound = LOW, HIGH
               rv(i, j, k) = rv(i, j, k) + outletFlag(ibound, idimn)*velforce(ibound, idimn)*outprofile(ibound, idimn)
            end do
         end do

         do ibound = LOW, HIGH
            QAuxLiq(ibound, JAXIS) = QAuxLiq(ibound, JAXIS) + &
                                     iLiq*outletFlag(ibound, JAXIS)*v(i, j, k)*volcell*outprofile(ibound, JAXIS)

            volAuxLiq(ibound, JAXIS) = volAuxLiq(ibound, JAXIS) + &
                                       iLiq*outletFlag(ibound, JAXIS)*volcell*outprofile(ibound, JAXIS)

            QAuxGas(ibound, JAXIS) = QAuxGas(ibound, JAXIS) + &
                                     iGas*outletFlag(ibound, JAXIS)*v(i, j, k)*volcell*outprofile(ibound, JAXIS)

            volAuxGas(ibound, JAXIS) = volAuxGas(ibound, JAXIS) + &
                                       iGas*outletFlag(ibound, JAXIS)*volcell*outprofile(ibound, JAXIS)
         end do

      end do
   end do

end subroutine sim_outletVelBlk2dPhased

subroutine sim_outletVelBlk3dPhased(u, v, w, ru, rv, rw, phi, xcenter, ycenter, zcenter, &
                                    boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                    inletFlag, inletBuffer, inletGrowthRate, &
                                    outletFlag, outletBuffer, outletGrowthRate, &
                                    volAuxLiq, volAuxGas, QAuxLiq, QAuxGas, QOutLiq, QOutGas, &
                                    xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ)

   implicit none
   real, dimension(:, :, :), intent(in)  :: u, v, w
   real, dimension(:, :, :), intent(inout) :: ru, rv, rw
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:), intent(in)        :: xcenter, ycenter, zcenter
   real, dimension(:, :), intent(in)     :: boundBox
   real, intent(in)                      :: dt, dx, dy, dz
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
   integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag, inletFlag
   real, intent(in) :: inletBuffer, inletGrowthRate, outletBuffer, outletGrowthRate
   real, intent(inout) :: QAuxLiq(LOW:HIGH, MDIM), QAuxGas(LOW:HIGH, MDIM)
   real, intent(inout) :: volAuxLiq(LOW:HIGH, MDIM), volAuxGas(LOW:HIGH, MDIM)
   real, intent(in) :: QOutLiq(LOW:HIGH, MDIM), QOutGas(LOW:HIGH, MDIM)
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ

   integer :: i, j, k

   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2 + 1
         end do
      end do
   end do

   do k = kz1, kz2
      do j = jy1, jy2 + 1
         do i = ix1, ix2
         end do
      end do
   end do

   do k = kz1, kz2 + 1
      do j = jy1, jy2
         do i = ix1, ix2
         end do
      end do
   end do

end subroutine sim_outletVelBlk3dPhased
