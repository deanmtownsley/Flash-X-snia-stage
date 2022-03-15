!!***if* source/Simulation/SimulationForcing/incompFlow/Outflow/sim_outflowVelBlk
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

subroutine sim_outflowVelBlk2d(velOut, u, v, ru, rv, xcell, ycell, &
                               boundBox, dt, dx, dy, ix1, ix2, jy1, jy2)

   use Simulation_data
   use sim_outflowData

   implicit none
   real, intent(inout)                     :: velOut
   real, dimension(:, :, :), intent(in)    :: u, v
   real, dimension(:, :, :), intent(inout) :: ru, rv
   real, dimension(:), intent(in)          :: xcell, ycell
   real, dimension(:, :), intent(in)       :: boundBox
   real, intent(in)                        :: dt, dx, dy
   integer, intent(in)                     :: ix1, ix2, jy1, jy2

   integer :: i, j, k
   real    :: xi, yi
   real    :: uyplus, uymins, vyplus, vymins

   k = 1

   do j = jy1 + 1, jy2 - 1
      do i = ix1 + 1, ix2
         xi = xcell(i)
         yi = ycell(j)

         uyplus = (u(i, j + 1, k) + u(i, j, k))*0.5
         uymins = (u(i, j, k) + u(i, j - 1, k))*0.5

         ru(i, j, k) = -sim_outflowVel(HIGH, JAXIS)*(uyplus - uymins)/dy* &
                       (1/(1 + exp(-sim_outflowGrowthRate*(yi - sim_outflowRegion))))
      end do
   end do

   do j = jy1 + 1, jy2
      do i = ix1 + 1, ix2 - 1
         xi = xcell(i)
         yi = ycell(j)

         vyplus = (v(i, j + 1, k) + v(i, j, k))*0.5
         vymins = (v(i, j, k) + v(i, j - 1, k))*0.5

         rv(i, j, k) = -sim_outflowVel(HIGH, JAXIS)*(vyplus - vymins)/dy* &
                       (1/(1 + exp(-sim_outflowGrowthRate*(yi - sim_outflowRegion))))

         if (yi .le. sim_yMax .and. yi .ge. sim_yMax - dy) velOut = max(velOut, v(i, j + 1, k))

      end do
   end do

end subroutine sim_outflowVelBlk2d

subroutine sim_outflowVelBlk3d(velOut, u, v, w, ru, rv, rw, xcell, ycell, zcell, &
                               boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)

   use Simulation_data
   use sim_outflowData

   implicit none
   real, intent(inout)                   :: velOut
   real, dimension(:, :, :), intent(in)  :: u, v, w
   real, dimension(:, :, :), intent(inout) :: ru, rv, rw
   real, dimension(:), intent(in)        :: xcell, ycell, zcell
   real, dimension(:, :), intent(in)     :: boundBox
   real, intent(in)                      :: dt, dx, dy, dz
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2

   integer :: i, j, k
   real    :: xi, yi, zi
   real    :: uyplus, uymins, vyplus, vymins, wyplus, wymins

   do k = kz1 + 1, kz2 - 1
      do j = jy1 + 1, jy2 - 1
         do i = ix1 + 1, ix2
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            uyplus = (u(i, j + 1, k) + u(i, j, k))*0.5
            uymins = (u(i, j, k) + u(i, j - 1, k))*0.5

            ru(i, j, k) = -sim_outflowVel(HIGH, JAXIS)*(uyplus - uymins)/dy* &
                          (1/(1 + exp(-sim_outflowGrowthRate*(yi - sim_outflowRegion))))

         end do
      end do
   end do

   do k = kz1 + 1, kz2 - 1
      do j = jy1 + 1, jy2
         do i = ix1 + 1, ix2 - 1
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            vyplus = (v(i, j + 1, k) + v(i, j, k))*0.5
            vymins = (v(i, j, k) + v(i, j - 1, k))*0.5

            rv(i, j, k) = -sim_outflowVel(HIGH, JAXIS)*(vyplus - vymins)/dy* &
                          (1/(1 + exp(-sim_outflowGrowthRate*(yi - sim_outflowRegion))))

            if (yi .le. sim_yMax .and. yi .ge. sim_yMax - dy) velOut = max(velOut, v(i, j + 1, k))

         end do
      end do
   end do

   do k = kz1 + 1, kz2
      do j = jy1 + 1, jy2 - 1
         do i = ix1 + 1, ix2 - 1
            xi = xcell(i)
            yi = ycell(j)
            zi = zcell(k)

            wyplus = (w(i, j + 1, k) + w(i, j, k))*0.5
            wymins = (w(i, j, k) + w(i, j - 1, k))*0.5

            rw(i, j, k) = -sim_outflowVel(HIGH, JAXIS)*(wyplus - wymins)/dy* &
                          (1/(1 + exp(-sim_outflowGrowthRate*(yi - sim_outflowRegion))))
         end do
      end do
   end do

end subroutine sim_outflowVelBlk3d
