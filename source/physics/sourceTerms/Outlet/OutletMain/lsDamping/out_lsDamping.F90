!!***if* source/physics/sourceTerms/Outlet/OutletMain/lsDamping/out_lsDamping
!!
!!
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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

subroutine out_lsDamping(pfrc, phi, xcenter, ycenter, zcenter, boundBox, &
                               dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                               outletFlag, outletSink, outletBuffer, outletGrowthRate, &
                               xMin, xMax, yMin, yMax, zMin, zMax)

   implicit none

   real, dimension(:, :, :), intent(inout) :: pfrc
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:), intent(in) :: xcenter, ycenter, zcenter
   real, dimension(:, :), intent(in) :: boundBox
   real, intent(in) :: dt, dx, dy, dz
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
   integer, dimension(2, MDIM), intent(in) :: outletFlag
   real, intent(in) :: outletSink, outletBuffer, outletGrowthRate
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax

   integer :: i, j, k
   real    :: xi, yi, zi

   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2
            xi = xcenter(i)
            yi = ycenter(j)
            zi = zcenter(k)

            pfrc(i, j, k) = pfrc(i, j, k) &
                            + outletSink*outletFlag(HIGH, IAXIS)*(2/(1 + exp(-outletGrowthRate*(xi - xMax)/outletBuffer))) &
                            + outletSink*outletFlag(HIGH, JAXIS)*(2/(1 + exp(-outletGrowthRate*(yi - yMax)/outletBuffer)))

#if NDIM == MDIM
            pfrc(i, j, k) = pfrc(i, j, k) &
                            + outletSink*outletFlag(HIGH, KAXIS)*(2/(1 + exp(-outletGrowthRate*(zi - zMax)/outletBuffer)))
#endif

         end do
      end do
   end do

end subroutine out_lsDamping
