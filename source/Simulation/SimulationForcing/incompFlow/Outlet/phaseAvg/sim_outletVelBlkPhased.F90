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
                                    boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, inletFlag, &
                                    outletFlag, phaseAuxLiq, phaseAuxGas, &
                                    velAuxLiq, velAuxGas, velOutLiq, velOutGas, &
                                    outletBuffer, outletGrowthRate, &
                                    xMin, xMax, yMin, yMax, gravX, gravY)

   implicit none
   real, dimension(:, :, :), intent(in)    :: u, v
   real, dimension(:, :, :), intent(inout) :: ru, rv
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:), intent(in)          :: xcell, ycell
   real, dimension(:, :), intent(in)       :: boundBox
   real, intent(in)                        :: dt, dx, dy
   integer, intent(in)                     :: ix1, ix2, jy1, jy2
   integer, dimension(2, MDIM), intent(in) :: outletFlag, inletFlag
   real, intent(inout) :: velAuxLiq(2, MDIM), velAuxGas(2, MDIM), phaseAuxLiq(2, MDIM), phaseAuxGas(2, MDIM)
   real, intent(in) :: velOutLiq(2, MDIM), velOutGas(2, MDIM), outletBuffer, outletGrowthRate
   real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY

   integer :: i, j, k
   real    :: xi, yi
   real    :: uplus, umins, vplus, vmins
   real    :: uforce, vforce
   real    :: phicell

   k = 1

   do j = jy1, jy2
      do i = ix1, ix2 + 1
         xi = xcell(i)
         yi = ycell(j)

         phicell = (phi(i, j, k) + phi(i - 1, j, k))*.5

         uforce = ((1 - sign(1., phicell))/2)*(velOutLiq(HIGH, IAXIS) - u(i, j, k))/dt &
                  + ((1 + sign(1., phicell))/2)*(velOutGas(HIGH, IAXIS) - u(i, j, k))/dt

         ru(i, j, k) = ru(i, j, k) &
                       + outletFlag(HIGH, IAXIS)*uforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                       + outletFlag(HIGH, JAXIS)*uforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer)))

         velAuxLiq(HIGH, IAXIS) = velAuxLiq(HIGH, IAXIS) &
                                  + outletFlag(HIGH, IAXIS)*u(i, j, k)*((1 - sign(1., phicell))/2) &
                                  *(dx/(xMax - xMin))*(dy/(yMax - yMin))

         phaseAuxLiq(HIGH, IAXIS) = phaseAuxLiq(HIGH, IAXIS) &
                                    + outletFlag(HIGH, IAXIS)*((1 - sign(1., phicell))/2) &
                                    *(dx/(xMax - xMin))*(dy/(yMax - yMin))

         velAuxGas(HIGH, IAXIS) = velAuxGas(HIGH, IAXIS) &
                                  + outletFlag(HIGH, IAXIS)*u(i, j, k)*((1 + sign(1., phicell))/2) &
                                  *(dx/(xMax - xMin))*(dy/(yMax - yMin))

         phaseAuxGas(HIGH, IAXIS) = phaseAuxGas(HIGH, IAXIS) &
                                    + outletFlag(HIGH, IAXIS)*((1 + sign(1., phicell))/2) &
                                    *(dx/(xMax - xMin))*(dy/(yMax - yMin))

      end do
   end do

   do j = jy1, jy2 + 1
      do i = ix1, ix2
         xi = xcell(i)
         yi = ycell(j)

         phicell = (phi(i, j, k) + phi(i, j - 1, k))*.5

         vforce = ((1 - sign(1., phicell))/2)*(velOutLiq(HIGH, JAXIS) - v(i, j, k))/dt &
                  + ((1 + sign(1., phicell))/2)*(velOutGas(HIGH, JAXIS) - v(i, j, k))/dt

         rv(i, j, k) = rv(i, j, k) &
                       + outletFlag(HIGH, IAXIS)*vforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                       + outletFlag(HIGH, JAXIS)*vforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer)))

         velAuxLiq(HIGH, JAXIS) = velAuxLiq(HIGH, JAXIS) &
                                  + outletFlag(HIGH, JAXIS)*v(i, j, k)*((1 - sign(1., phicell))/2) &
                                  *(dx/(xMax - xMin))*(dy/(yMax - yMin))

         phaseAuxLiq(HIGH, JAXIS) = phaseAuxLiq(HIGH, JAXIS) &
                                    + outletFlag(HIGH, JAXIS)*((1 - sign(1., phicell))/2) &
                                    *(dx/(xMax - xMin))*(dy/(yMax - yMin))

         velAuxGas(HIGH, JAXIS) = velAuxGas(HIGH, JAXIS) &
                                  + outletFlag(HIGH, JAXIS)*v(i, j, k)*((1 + sign(1., phicell))/2) &
                                  *(dx/(xMax - xMin))*(dy/(yMax - yMin))

         phaseAuxGas(HIGH, JAXIS) = phaseAuxGas(HIGH, JAXIS) &
                                    + outletFlag(HIGH, JAXIS)*((1 + sign(1., phicell))/2) &
                                    *(dx/(xMax - xMin))*(dy/(yMax - yMin))

      end do
   end do

end subroutine sim_outletVelBlk2dPhased

subroutine sim_outletVelBlk3dPhased(u, v, w, ru, rv, rw, phi, xcell, ycell, zcell, &
                                    boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, inletFlag, &
                                    outletFlag, phaseAuxLiq, phaseAuxGas, &
                                    velAuxLiq, velAuxGas, velOutLiq, velOutGas, &
                                    outletBuffer, outletGrowthRate, &
                                    xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ)

   implicit none
   real, dimension(:, :, :), intent(in)  :: u, v, w
   real, dimension(:, :, :), intent(inout) :: ru, rv, rw
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:), intent(in)        :: xcell, ycell, zcell
   real, dimension(:, :), intent(in)     :: boundBox
   real, intent(in)                      :: dt, dx, dy, dz
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
   integer, dimension(2, MDIM), intent(in) :: outletFlag, inletFlag
   real, intent(inout) :: velAuxLiq(2, MDIM), velAuxGas(2, MDIM), phaseAuxLiq(2, MDIM), phaseAuxGas(2, MDIM)
   real, intent(in) :: velOutLiq(2, MDIM), velOutGas(2, MDIM), outletBuffer, outletGrowthRate
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ

   integer :: i, j, k
   real    :: xi, yi, zi
   real    :: uplus, umins, vplus, vmins, wplus, wmins
   real    :: uforce, vforce, wforce
   real    :: phicell

   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2 + 1
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            phicell = (phi(i, j, k) + phi(i - 1, j, k))*.5

            uforce = ((1 - sign(1., phicell))/2)*(velOutLiq(HIGH, IAXIS) - u(i, j, k))/dt &
                     + ((1 + sign(1., phicell))/2)*(velOutGas(HIGH, IAXIS) - u(i, j, k))/dt

            ru(i, j, k) = ru(i, j, k) &
                          + outletFlag(HIGH, IAXIS)*uforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                          + outletFlag(HIGH, JAXIS)*uforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer))) &
                          + outletFlag(HIGH, KAXIS)*uforce*(2/(1 + exp(-outletGrowthRate*(zi - zMax)/outletBuffer)))

            velAuxLiq(HIGH, IAXIS) = velAuxLiq(HIGH, IAXIS) &
                                     + outletFlag(HIGH, IAXIS)*u(i, j, k)*((1 - sign(1., phicell))/2) &
                                     *(dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

            phaseAuxLiq(HIGH, IAXIS) = phaseAuxLiq(HIGH, IAXIS) &
                                       + outletFlag(HIGH, IAXIS)*((1 - sign(1., phicell))/2) &
                                       *(dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

            velAuxGas(HIGH, IAXIS) = velAuxGas(HIGH, IAXIS) &
                                     + outletFlag(HIGH, IAXIS)*u(i, j, k)*((1 + sign(1., phicell))/2) &
                                     *(dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

            phaseAuxGas(HIGH, IAXIS) = phaseAuxGas(HIGH, IAXIS) &
                                       + outletFlag(HIGH, IAXIS)*((1 + sign(1., phicell))/2) &
                                       *(dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

         end do
      end do
   end do

   do k = kz1, kz2
      do j = jy1, jy2 + 1
         do i = ix1, ix2
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            phicell = (phi(i, j, k) + phi(i, j - 1, k))*.5

            vforce = ((1 - sign(1., phicell))/2)*(velOutLiq(HIGH, JAXIS) - v(i, j, k))/dt &
                     + ((1 + sign(1., phicell))/2)*(velOutGas(HIGH, JAXIS) - v(i, j, k))/dt

            rv(i, j, k) = rv(i, j, k) &
                          + outletFlag(HIGH, IAXIS)*vforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                          + outletFlag(HIGH, JAXIS)*vforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer))) &
                          + outletFlag(HIGH, KAXIS)*vforce*(2/(1 + exp(-outletGrowthRate*(zi - zMax)/outletBuffer)))

            velAuxLiq(HIGH, JAXIS) = velAuxLiq(HIGH, JAXIS) &
                                     + outletFlag(HIGH, JAXIS)*v(i, j, k)*((1 - sign(1., phicell))/2) &
                                     *(dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

            phaseAuxLiq(HIGH, JAXIS) = phaseAuxLiq(HIGH, JAXIS) &
                                       + outletFlag(HIGH, JAXIS)*((1 - sign(1., phicell))/2) &
                                       *(dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

            velAuxGas(HIGH, JAXIS) = velAuxGas(HIGH, JAXIS) &
                                     + outletFlag(HIGH, JAXIS)*v(i, j, k)*((1 + sign(1., phicell))/2) &
                                     *(dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

            phaseAuxGas(HIGH, JAXIS) = phaseAuxGas(HIGH, JAXIS) &
                                       + outletFlag(HIGH, JAXIS)*((1 + sign(1., phicell))/2) &
                                       *(dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

         end do
      end do
   end do

   do k = kz1, kz2 + 1
      do j = jy1, jy2
         do i = ix1, ix2
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            phicell = (phi(i, j, k) + phi(i, j, k - 1))*.5

            wforce = ((1 - sign(1., phicell))/2)*(velOutLiq(HIGH, KAXIS) - w(i, j, k))/dt &
                     + ((1 + sign(1., phicell))/2)*(velOutGas(HIGH, KAXIS) - w(i, j, k))/dt

            rw(i, j, k) = rw(i, j, k) &
                          + outletFlag(HIGH, IAXIS)*wforce*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                          + outletFlag(HIGH, JAXIS)*wforce*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer))) &
                          + outletFlag(HIGH, KAXIS)*wforce*(2/(1 + exp(-outletGrowthRate*(zi - zMax)/outletBuffer)))

            velAuxLiq(HIGH, KAXIS) = velAuxLiq(HIGH, KAXIS) &
                                     + outletFlag(HIGH, KAXIS)*w(i, j, k)*((1 - sign(1., phicell))/2) &
                                     *(dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

            phaseAuxLiq(HIGH, KAXIS) = phaseAuxLiq(HIGH, KAXIS) &
                                       + outletFlag(HIGH, KAXIS)*((1 - sign(1., phicell))/2) &
                                       *(dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

            velAuxGas(HIGH, KAXIS) = velAuxGas(HIGH, KAXIS) &
                                     + outletFlag(HIGH, KAXIS)*w(i, j, k)*((1 + sign(1., phicell))/2) &
                                     *(dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

            phaseAuxGas(HIGH, KAXIS) = phaseAuxGas(HIGH, KAXIS) &
                                       + outletFlag(HIGH, KAXIS)*((1 + sign(1., phicell))/2) &
                                       *(dx/(xMax - xMin))*(dy/(yMax - yMin))*(dz/(zMax - zMin))

         end do
      end do
   end do

end subroutine sim_outletVelBlk3dPhased
