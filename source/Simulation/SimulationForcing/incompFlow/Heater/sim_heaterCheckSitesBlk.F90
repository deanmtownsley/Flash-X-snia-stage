!!***if* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterCheckSitesBlk
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

subroutine sim_heaterCheckSitesBlk2d(phi, xcell, ycell, boundBox, ix1, ix2, jy1, jy2)

   use sim_heaterData

   implicit none
   real, dimension(:, :, :), intent(in)      :: phi
   real, dimension(:), intent(in)          :: xcell, ycell
   real, dimension(:, :), intent(in)        :: boundBox
   integer, intent(in)                    :: ix1, ix2, jy1, jy2

   integer :: i, j, k, isite, htr
   real    :: xi, xp, yi, yp
   real    :: phiSW, phiSE, phiNW, phiNE, phiSite

   type(sim_heaterType), pointer :: heater

   k = 1

   do htr = 1, sim_numHeaters

      heater => sim_heaterInfo(htr)

      if (boundBox(HIGH, IAXIS) .le. heater%xMin .or. boundBox(LOW, IAXIS) .ge. heater%xMax .or. &
          boundBox(HIGH, JAXIS) .le. heater%yMin .or. boundBox(LOW, JAXIS) .ge. heater%yMax) cycle

      do j = jy1, jy2 - 1
         do i = ix1, ix2 - 1
            do isite = 1, heater%numSites
               xi = xcell(i)
               xp = xcell(i + 1)
               yi = ycell(j)
               yp = ycell(j + 1)

               phiSW = phi(i, j, k)
               phiSE = phi(i + 1, j, k)
               phiNW = phi(i, j + 1, k)
               phiNE = phi(i + 1, j + 1, k)

               phiSite = (phiSW + phiSE + phiNW + phiNE)/4.

               if (xi .le. heater%xSite(isite) .and. xp .ge. heater%xSite(isite) .and. &
                   yi .le. heater%ySite(isite) .and. yp .ge. heater%ySite(isite)) then

                  if (phiSite .ge. 0.0) then
                     heater%siteIsAttachedCurr(isite) = heater%siteIsAttachedCurr(isite) .or. .true.
                  else
                     heater%siteIsAttachedCurr(isite) = heater%siteIsAttachedCurr(isite) .or. .false.
                  end if

               end if

            end do
         end do
      end do

   end do

   return
end subroutine sim_heaterCheckSitesBlk2d

subroutine sim_heaterCheckSitesBlk3d(phi, xcell, ycell, zcell, boundBox, ix1, ix2, jy1, jy2, kz1, kz2)

   use sim_heaterData

   implicit none
   real, dimension(:, :, :), intent(in)  :: phi
   real, dimension(:), intent(in)      :: xcell, ycell, zcell
   real, dimension(:, :), intent(in)    :: boundBox
   integer, intent(in)                :: ix1, ix2, jy1, jy2, kz1, kz2

   integer :: i, j, k, isite, htr
   real    :: xi, xp, yi, yp, zi, zp
   real    :: phiFSW, phiFSE, phiFNW, phiFNE
   real    :: phiBSW, phiBSE, phiBNW, phiBNE
   real    :: phiSite

   type(sim_heaterType), pointer :: heater

   do htr = 1, sim_numHeaters

      heater => sim_heaterInfo(htr)

      if (boundBox(HIGH, IAXIS) .le. heater%xMin .or. boundBox(LOW, IAXIS) .ge. heater%xMax .or. &
          boundBox(HIGH, JAXIS) .le. heater%yMin .or. boundBox(LOW, JAXIS) .ge. heater%yMax .or. &
          boundBox(HIGH, KAXIS) .le. heater%zMin .or. boundBox(LOW, KAXIS) .ge. heater%zMax) cycle

      do k = kz1, kz2 - 1
         do j = jy1, jy2 - 1
            do i = ix1, ix2 - 1
               do isite = 1, heater%numSites

                  xi = xcell(i)
                  xp = xcell(i + 1)
                  yi = ycell(j)
                  yp = ycell(j + 1)
                  zi = zcell(k)
                  zp = zcell(k + 1)

                  phiFSW = phi(i, j, k)
                  phiFSE = phi(i + 1, j, k)
                  phiFNW = phi(i, j + 1, k)
                  phiFNE = phi(i + 1, j + 1, k)

                  phiBSW = phi(i, j, k + 1)
                  phiBSE = phi(i + 1, j, k + 1)
                  phiBNW = phi(i, j + 1, k + 1)
                  phiBNE = phi(i + 1, j + 1, k + 1)

                  phiSite = (phiFSW + phiFSE + phiFNW + phiFNE + phiBSW + phiBSE + phiBNW + phiBNE)/8.

                  if (xi .le. heater%xSite(isite) .and. xp .ge. heater%xSite(isite) .and. &
                      yi .le. heater%ySite(isite) .and. yp .ge. heater%ySite(isite) .and. &
                      zi .le. heater%zSite(isite) .and. zp .ge. heater%zSite(isite)) then

                     if (phiSite .ge. 0.0) then
                        heater%siteIsAttachedCurr(isite) = heater%siteIsAttachedCurr(isite) .or. .true.
                     else
                        heater%siteIsAttachedCurr(isite) = heater%siteIsAttachedCurr(isite) .or. .false.
                     end if

                  end if

               end do
            end do
         end do
      end do

   end do

   return
end subroutine sim_heaterCheckSitesBlk3d
