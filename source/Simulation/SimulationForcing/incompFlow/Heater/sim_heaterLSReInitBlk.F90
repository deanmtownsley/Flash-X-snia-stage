!!***if* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterLSReInitBlk
!!
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

subroutine sim_heaterLSReInitBlk(phi, xcell, ycell, zcell, boundBox, stime, ix1, ix2, jy1, jy2, kz1, kz2)

   use sim_heaterData

   implicit none
   real, dimension(:, :, :), intent(inout)  :: phi
   real, dimension(:), intent(in)        :: xcell, ycell, zcell
   real, dimension(:, :), intent(in)      :: boundBox
   real, intent(in)                      :: stime
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2

   type(sim_heaterType), pointer  :: heater
   integer :: i, j, k, htr, isite
   real    :: idfun, iseedY, iseedX, iseedZ, iradius

   do htr = 1, sim_numHeaters

      heater => sim_heaterInfo(htr)

      if (boundBox(HIGH, IAXIS) .le. heater%xMin .or. boundBox(LOW, IAXIS) .ge. heater%xMax .or. &
          boundBox(HIGH, JAXIS) .le. heater%yMin .or. boundBox(LOW, JAXIS) .ge. heater%yMax .or. &
          boundBox(HIGH, KAXIS) .le. heater%zMin .or. boundBox(LOW, KAXIS) .ge. heater%zMax) cycle

      do k = kz1, kz2
         do j = jy1, jy2
            do i = ix1, ix2
               do isite = 1, heater%numSites

                  if (((heater%siteTimeStamp(isite) + heater%nucWaitTime) .le. stime) .and. &
                      (heater%siteIsAttachedPrev(isite) .eqv. .false.)) then
                     iradius = heater%seedRadius
                     iseedX = heater%xSite(isite)
                     iseedZ = heater%zSite(isite)
                     iseedY = heater%ySite(isite) + heater%seedHeight
                     idfun = iradius - sqrt((xcell(i) - iseedX)**2 + (ycell(j) - iseedY)**2 + (zcell(k) - iseedZ)**2)
                     phi(i, j, k) = max(phi(i, j, k), idfun)
                  end if

               end do
            end do
         end do
      end do

   end do

   return

end subroutine sim_heaterLSReInitBlk
