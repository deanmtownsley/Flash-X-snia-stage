!!***if* source/Simulation/SimulationForcing/incompFlow/Outlet/sim_outletVelFrc
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

subroutine sim_outletVelFrc(vel, rhs, xgrid, ygrid, zgrid, &
                            dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                            xMin, xMax, yMin, yMax, zMin, zMax, &
                            outletFlag, outletBuffer, outletGrowthRate, &
                            axis, volAux, QAux, QOut)

   implicit none
   real, dimension(:, :, :), intent(in) :: vel
   real, dimension(:, :, :), intent(inout) :: rhs
   real, dimension(:), intent(in) :: xgrid, ygrid, zgrid
   real, intent(in) :: dt, dx, dy, dz
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax
   integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag
   real, intent(in) :: outletBuffer, outletGrowthRate
   integer, intent(in) :: axis
   real, intent(inout) :: QAux(LOW:HIGH, MDIM), volAux(LOW:HIGH, MDIM)
   real, intent(in) :: QOut(LOW:HIGH, MDIM)

   !---Local variables
   integer :: i, j, k, idimn, ibound, iforce, inorm
   real :: xcell, ycell, zcell, velout, velforce
   real :: outprofile(LOW:HIGH, MDIM), velgrad(LOW:HIGH, MDIM)
   real, parameter :: velref = 1.

   ! Loop over indices
   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2

            ! Get cell co-ordinates from grid array
            xcell = xgrid(i)
            ycell = ygrid(j)
            zcell = zgrid(k)

            ! Setup outflow profile using a sigmod function
            outprofile(LOW, IAXIS:JAXIS) = (/2/(1 + exp(outletGrowthRate*(xcell - xMin)/outletBuffer)), &
                                             2/(1 + exp(outletGrowthRate*(ycell - yMin)/outletBuffer))/)
            !
            outprofile(HIGH, IAXIS:JAXIS) = (/2/(1 + exp(-outletGrowthRate*(xcell - xMax)/outletBuffer)), &
                                              2/(1 + exp(-outletGrowthRate*(ycell - yMax)/outletBuffer))/)

            ! Calculate gradient normal to the boundary
            velgrad(LOW, IAXIS:JAXIS) = (/(vel(i - 1, j, k) - vel(i + 1, j, k))/(2*dx), &
                                          (vel(i, j - 1, k) - vel(i, j + 1, k))/(2*dy)/)
            !
            velgrad(HIGH, IAXIS:JAXIS) = (/(vel(i + 1, j, k) - vel(i - 1, j, k))/(2*dx), &
                                           (vel(i, j + 1, k) - vel(i, j - 1, k))/(2*dy)/)

#if NDIM == MDIM
            outprofile(LOW, KAXIS) = 2/(1 + exp(outletGrowthRate*(zcell - zMin)/outletBuffer))
            outprofile(HIGH, KAXIS) = 2/(1 + exp(-outletGrowthRate*(zcell - zMax)/outletBuffer))

            velgrad(LOW, KAXIS) = (vel(i, j, k - 1) - vel(i, j, k + 1))/(2*dz)
            velgrad(HIGH, KAXIS) = (vel(i, j, k + 1) - vel(i, j, k - 1))/(2*dz)
#endif

            ! Loop over all boundaries
            do idimn = 1, NDIM
               do ibound = LOW, HIGH

                  ! Get absolute outlet velocity
                  ! QOut is the mean outlet velocity
                  ! velref is equal to 1 (reference scale)
                  velout = max(velref, abs(QOut(ibound, idimn)))

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

            ! Update QAux and volAux on local processor
            ! This is later used to compute QOut in Simulation_adjustEvolution
            do ibound = LOW, HIGH
               QAux(ibound, axis) = QAux(ibound, axis) + outletFlag(ibound, axis)*vel(i, j, k)*outprofile(ibound, axis)
               volAux(ibound, axis) = volAux(ibound, axis) + outletFlag(ibound, axis)*outprofile(ibound, axis)
            end do

         end do
      end do
   end do

end subroutine sim_outletVelFrc
