!!***if* source/Simulation/SimulationMain/incompFlow/CounterFlow/sim_inletVelBlkPhased
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

subroutine sim_inletVelBlk2dPhased(u, v, ru, rv, phi, xcell, ycell, &
                                   boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, &
                                   inletFlag, outletFlag, &
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
   real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY

   integer :: i, j, k
   real :: xi, yi, phicell
   real :: uforce, vforce
   real, parameter :: inletGrowthRate = 4.0, inletBuffer = 1.0

   k = 1

   do j = jy1, jy2
      do i = ix1, ix2 + 1
         xi = xcell(i)
         yi = ycell(j)

         phicell = (phi(i, j, k) + phi(i - 1, j, k))*.5

         uforce = ((1 - sign(1., phicell))/2)*(-1.0 - u(i, j, k))/dt &
                  + ((1 + sign(1., phicell))/2)*(1.0 - u(i, j, k))/dt

         ru(i, j, k) = ru(i, j, k) + uforce*(2/(1 + exp(inletGrowthRate*(xi - xMin)/inletBuffer)))

      end do
   end do

   do j = jy1, jy2 + 1
      do i = ix1, ix2
         xi = xcell(i)
         yi = ycell(j)

         phicell = (phi(i, j, k) + phi(i, j - 1, k))*.5

         vforce = ((1 - sign(1., phicell))/2)*(0.-v(i, j, k))/dt &
                  + ((1 + sign(1., phicell))/2)*(0.-v(i, j, k))/dt

         rv(i, j, k) = rv(i, j, k) + vforce*(2/(1 + exp(inletGrowthRate*(xi - xMin)/inletBuffer)))

      end do
   end do

end subroutine sim_inletVelBlk2dPhased

subroutine sim_inletVelBlk3dPhased(u, v, w, ru, rv, rw, phi, xcell, ycell, zcell, &
                                   boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                   inletFlag, outletFlag, &
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
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ

end subroutine sim_inletVelBlk3dPhased
