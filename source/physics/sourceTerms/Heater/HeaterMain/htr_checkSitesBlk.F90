!!***if* source/physics/sourceTerms/Heater/HeaterMain/htr_checkSitesBlk
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

subroutine htr_checkSitesBlk2d(phi, xcell, ycell, boundBox, ix1, ix2, jy1, jy2, lblock)

   use Heater_data
   use Heater_type, ONLY: Heater_type_t
   use Driver_interface, ONLY: Driver_abort

   implicit none
   real, dimension(:, :, :), intent(in)      :: phi
   real, dimension(:), intent(in)          :: xcell, ycell
   real, dimension(:, :), intent(in)        :: boundBox
   integer, intent(in)                    :: ix1, ix2, jy1, jy2, lblock

   integer :: i, j, k, isite, htr, annIndex, isiteblk
   real    :: xi, xp, yi, yp
   real    :: phiSW, phiSE, phiNW, phiNE, phiSite

   type(Heater_type_t), pointer :: heater

   k = 1

   do htr = 1, htr_numHeaters

      heater => htr_heaterInfo(htr)

      if (boundBox(HIGH, IAXIS) .le. heater%xMin .or. boundBox(LOW, IAXIS) .ge. heater%xMax .or. &
          boundBox(HIGH, JAXIS) .le. heater%yMin .or. boundBox(LOW, JAXIS) .ge. heater%yMax) cycle

      do j = jy1, jy2-1
         do i = ix1, ix2-1
            do isiteblk = 1, heater%numSitesBlk(lblock)
               isite = heater%siteMapOnProc(lblock, isiteblk)

               if (isite < 1) call Driver_abort("[htr_checkSitesBlk2d] isite < 1")

               xi = xcell(i)
               xp = xcell(i+1)
               yi = ycell(j)
               yp = ycell(j+1)

               phiSW = phi(i, j, k)
               phiSE = phi(i+1, j, k)
               phiNW = phi(i, j+1, k)
               phiNE = phi(i+1, j+1, k)

               phiSite = (phiSW+phiSE+phiNW+phiNE)/4.

               if (xi .le. heater%xSiteProc(isite) .and. xp .ge. heater%xSiteProc(isite) .and. &
                   yi .le. heater%ySiteProc(isite) .and. yp .ge. heater%ySiteProc(isite)) then

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
end subroutine htr_checkSitesBlk2d

subroutine htr_checkSitesBlk3d(phi, xcell, ycell, zcell, boundBox, ix1, ix2, jy1, jy2, kz1, kz2, lblock)

   use Heater_data
   use Heater_type, ONLY: Heater_type_t
   use Driver_interface, ONLY: Driver_abort

   implicit none
   real, dimension(:, :, :), intent(in)  :: phi
   real, dimension(:), intent(in)      :: xcell, ycell, zcell
   real, dimension(:, :), intent(in)    :: boundBox
   integer, intent(in)                :: ix1, ix2, jy1, jy2, kz1, kz2, lblock

   integer :: i, j, k, isite, htr, isiteblk
   real    :: xi, xp, yi, yp, zi, zp
   real    :: phiFSW, phiFSE, phiFNW, phiFNE
   real    :: phiBSW, phiBSE, phiBNW, phiBNE
   real    :: phiSite

   type(Heater_type_t), pointer :: heater

   do htr = 1, htr_numHeaters

      heater => htr_heaterInfo(htr)

      if (boundBox(HIGH, IAXIS) .le. heater%xMin .or. boundBox(LOW, IAXIS) .ge. heater%xMax .or. &
          boundBox(HIGH, JAXIS) .le. heater%yMin .or. boundBox(LOW, JAXIS) .ge. heater%yMax .or. &
          boundBox(HIGH, KAXIS) .le. heater%zMin .or. boundBox(LOW, KAXIS) .ge. heater%zMax) cycle

      do k = kz1, kz2-1
         do j = jy1, jy2-1
            do i = ix1, ix2-1
               do isiteblk = 1, heater%numSitesBlk(lblock)
                  isite = heater%siteMapOnProc(lblock, isiteblk)

                  if (isite < 1) call Driver_abort("[htr_checkSitesBlk3d] isite < 1")

                  xi = xcell(i)
                  xp = xcell(i+1)
                  yi = ycell(j)
                  yp = ycell(j+1)
                  zi = zcell(k)
                  zp = zcell(k+1)

                  phiFSW = phi(i, j, k)
                  phiFSE = phi(i+1, j, k)
                  phiFNW = phi(i, j+1, k)
                  phiFNE = phi(i+1, j+1, k)

                  phiBSW = phi(i, j, k+1)
                  phiBSE = phi(i+1, j, k+1)
                  phiBNW = phi(i, j+1, k+1)
                  phiBNE = phi(i+1, j+1, k+1)

                  phiSite = (phiFSW+phiFSE+phiFNW+phiFNE+phiBSW+phiBSE+phiBNW+phiBNE)/8.

                  if (xi .le. heater%xSiteProc(isite) .and. xp .ge. heater%xSiteProc(isite) .and. &
                      yi .le. heater%ySiteProc(isite) .and. yp .ge. heater%ySiteProc(isite) .and. &
                      zi .le. heater%zSiteProc(isite) .and. zp .ge. heater%zSiteProc(isite)) then

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
end subroutine htr_checkSitesBlk3d
