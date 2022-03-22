!!***if* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterApplyBCToBlk
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

subroutine sim_heaterApplyBCToBlk2d(pfrc, tfrc, phi, temp, xcell, ycell, boundBox, dt, dx, dy, ix1, ix2, jy1, jy2)

   use Simulation_data
   use sim_heaterData

   implicit none
   real, dimension(:, :, :), intent(inout)   :: pfrc, tfrc
   real, dimension(:, :, :), intent(in)      :: phi, temp
   real, dimension(:), intent(in)          :: xcell, ycell
   real, dimension(:, :), intent(in)        :: boundBox
   real, intent(in)                        :: dt, dx, dy
   integer, intent(in)                    :: ix1, ix2, jy1, jy2

   integer :: i, j, k, htr
   real    :: xi, yi
   real    :: dynamicAngle, phiWall
   type(sim_heaterType), pointer :: heater

   k = 1

   do j = jy1, jy2
      do i = ix1, ix2

         xi = xcell(i)
         yi = ycell(j)

         do htr = 1, sim_numHeaters
            heater => sim_heaterInfo(htr)

         end do
      end do
   end do

   return
end subroutine sim_heaterApplyBCToBlk2d
