!!***if* source/Simulation/SimulationForcing/incompFlow/Outlet/sim_outletVelFrcPhased
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

subroutine sim_outletVelFrcPhased(vel, rhs, phi, xgrid, ygrid, zgrid, &
                                  dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                  xMin, xMax, yMin, yMax, zMin, zMax, &
                                  outletFlag, outletBuffer, outletGrowthRate, &
                                  axis, volAuxLiq, volAuxGas, QAuxLiq, QAuxGas, &
                                  QOutLiq, QOutGas)

   implicit none
   real, dimension(:, :, :), intent(in) :: vel, phi
   real, dimension(:, :, :), intent(inout) :: rhs
   real, dimension(:), intent(in) :: xgrid, ygrid, zgrid
   real, intent(in) :: dt, dx, dy, dz
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax
   integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag
   real, intent(in) :: outletBuffer, outletGrowthRate
   integer, intent(in) :: axis
   real, intent(inout) :: QAuxLiq(LOW:HIGH, MDIM), QAuxGas(LOW:HIGH, MDIM)
   real, intent(inout) :: volAuxLiq(LOW:HIGH, MDIM), volAuxGas(LOW:HIGH, MDIM)
   real, intent(in) :: QOutLiq(LOW:HIGH, MDIM), QOutGas(LOW:HIGH, MDIM)

   !---Local variables
   integer :: i, j, k, idimn, ibound, iforce, iliq, igas, inorm
   real :: xcell, ycell, zcell, velout, velforce
   real :: outprofile(LOW:HIGH, MDIM), velgrad(LOW:HIGH, MDIM), phiface(MDIM)
   real, parameter :: velref = 1.

   ! Loop over indices
   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2

            ! Get cell coordinates from grid arrays
            xcell = xgrid(i)
            ycell = ygrid(j)
            zcell = zgrid(k)

            ! Compute level set function on the cell face
            phiface(IAXIS:JAXIS) = (/(phi(i, j, k) + phi(i - 1, j, k))*.5, &
                                     (phi(i, j, k) + phi(i, j - 1, k))*.5/)

            ! Set outlet profile
            outprofile(LOW, IAXIS:JAXIS) = (/2/(1 + exp(outletGrowthRate*(xcell - xMin)/outletBuffer)), &
                                             2/(1 + exp(outletGrowthRate*(ycell - yMin)/outletBuffer))/)
            !
            outprofile(HIGH, IAXIS:JAXIS) = (/2/(1 + exp(-outletGrowthRate*(xcell - xMax)/outletBuffer)), &
                                              2/(1 + exp(-outletGrowthRate*(ycell - yMax)/outletBuffer))/)

            velgrad(LOW, IAXIS:JAXIS) = (/(vel(i - 1, j, k) - vel(i + 1, j, k))/(2*dx), &
                                          (vel(i, j - 1, k) - vel(i, j + 1, k))/(2*dy)/)
            !
            velgrad(HIGH, IAXIS:JAXIS) = (/(vel(i + 1, j, k) - vel(i - 1, j, k))/(2*dx), &
                                           (vel(i, j + 1, k) - vel(i, j - 1, k))/(2*dy)/)

#if NDIM == MDIM
            phiface(:) = (phi(i, j, k) + phi(i, j, k - 1))*.5

            outprofile(LOW, KAXIS) = 2/(1 + exp(outletGrowthRate*(zcell - zMin)/outletBuffer))
            outprofile(HIGH, KAXIS) = 2/(1 + exp(-outletGrowthRate*(zcell - zMax)/outletBuffer))

            velgrad(LOW, KAXIS) = (vel(i, j, k - 1) - vel(i, j, k + 1))/(2*dz)
            velgrad(HIGH, KAXIS) = (vel(i, j, k + 1) - vel(i, j, k - 1))/(2*dz)
#endif

            ! Loop over boundaries
            do idimn = 1, NDIM
               do ibound = LOW, HIGH

                  ! Cache the fluid information liquid or gas
                  iliq = (1 - int(sign(1., phiface(idimn))))/2
                  igas = (1 + int(sign(1., phiface(idimn))))/2

                  ! Get absolute outlet velocity
                  ! QOut is the mean outlet velocity
                  ! velref is equal to 1 (reference scale)
                  velout = max(velref, iliq*abs(QOutLiq(ibound, idimn)) + igas*abs(QOutGas(ibound, idimn)))

                  ! Check if local velocity greater than
                  ! outlet velocity
                  iforce = 0
                  if (abs(vel(i, j, k)) > velout) iforce = 1

                  ! Check if normal axis
                  inorm = 0
                  if (axis == idimn) inorm = 1

                  ! compute forcing on the local cell
                  velforce = inorm*iforce*((velout*vel(i, j, k)/(abs(vel(i, j, k)) + 1e-13) - vel(i, j, k))/dt - &
                                           velref*velgrad(ibound, idimn)) - (1 - inorm)*vel(i, j, k)/dt

                  ! Set source term for navier-stokes equation
                  rhs(i, j, k) = rhs(i, j, k) + velforce*outletFlag(ibound, idimn)*outprofile(ibound, idimn)
               end do
            end do

            ! Set QAux and volAux for liquid and gas phase
            ! will be used in Simulation_adjustEvolution to compute
            ! QOut
            do ibound = LOW, HIGH
               iliq = (1 - int(sign(1., phiface(axis))))/2
               igas = (1 + int(sign(1., phiface(axis))))/2

               QAuxLiq(ibound, axis) = QAuxLiq(ibound, axis) + iliq*outletFlag(ibound, axis)*vel(i, j, k)*outprofile(ibound, axis)
               volAuxLiq(ibound, axis) = volAuxLiq(ibound, axis) + iliq*outletFlag(ibound, axis)*outprofile(ibound, axis)

               QAuxGas(ibound, axis) = QAuxGas(ibound, axis) + igas*outletFlag(ibound, axis)*vel(i, j, k)*outprofile(ibound, axis)
               volAuxGas(ibound, axis) = volAuxGas(ibound, axis) + igas*outletFlag(ibound, axis)*outprofile(ibound, axis)
            end do

         end do
      end do
   end do

end subroutine sim_outletVelFrcPhased
