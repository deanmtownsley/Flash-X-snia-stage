!!***if* source/physics/sourceTerms/Heater/HeaterMain/htr_lsReInitBlk
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

subroutine htr_lsReInitBlk(phi, xcell, ycell, zcell, boundBox, stime, ix1, ix2, jy1, jy2, kz1, kz2, lblock)

   use Simulation_data, ONLY: sim_yMin, sim_yMax     
   use Heater_data
   use Heater_type, ONLY: Heater_type_t
   use Driver_interface, ONLY: Driver_abort
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get

   implicit none
   real, dimension(:, :, :), intent(inout)  :: phi
   real, dimension(:), intent(in)        :: xcell, ycell, zcell
   real, dimension(:, :), intent(in)      :: boundBox
   real, intent(in)                      :: stime
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2, lblock
   type(Heater_type_t), pointer  :: heater
   integer :: i, j, k, htr, isite, annIndex, isiteblk
   real    :: idfun, iseedY, iseedX, iseedZ, iradius

   call RuntimeParameters_get('ymin', sim_yMin)
   call RuntimeParameters_get('ymax', sim_yMax)

   do htr = 1, htr_numHeaters

      heater => htr_heaterInfo(htr)

#if NDIM < MDIM
      if (boundBox(HIGH, IAXIS) .le. heater%xMin .or. boundBox(LOW, IAXIS) .ge. heater%xMax .or. &
          boundBox(HIGH, JAXIS) .le. heater%yMin .or. boundBox(LOW, JAXIS) .ge. heater%yMax) cycle
#else
      if (boundBox(HIGH, IAXIS) .le. heater%xMin .or. boundBox(LOW, IAXIS) .ge. heater%xMax .or. &
          boundBox(HIGH, JAXIS) .le. heater%yMin .or. boundBox(LOW, JAXIS) .ge. heater%yMax .or. &
          boundBox(HIGH, KAXIS) .le. heater%zMin .or. boundBox(LOW, KAXIS) .ge. heater%zMax) cycle
#endif

      do k = kz1, kz2
         do j = jy1, jy2
            do i = ix1, ix2
               do isiteblk = 1, heater%numSitesBlk(lblock)

                  isite = heater%siteMapOnProc(lblock, isiteblk)
                  if (isite < 1) call Driver_abort("[htr_lsReInitBlk] isite < 1")

                  if (((heater%siteTimeStamp(isite)+heater%nucWaitTime) .le. stime) .and. &
                      (heater%siteIsAttachedPrev(isite) .eqv. .false.)) then
                     iradius = heater%seedRadius
                     iseedX = heater%xSiteProc(isite)
                     iseedZ = heater%zSiteProc(isite)
                     if(  abs(heater%ySiteProc(isite) - sim_yMin) .lt. abs(heater%ySiteProc(isite) - sim_yMax))then
                        iseedY = heater%ySiteProc(isite)+heater%seedHeight
                     else
                        iseedY = heater%ySiteProc(isite)-heater%seedHeight                             
                     end if
                     idfun = iradius-sqrt((xcell(i)-iseedX)**2+(ycell(j)-iseedY)**2+(zcell(k)-iseedZ)**2)
                     phi(i, j, k) = max(phi(i, j, k), idfun)
                  end if

               end do
            end do
         end do
      end do
   end do
end subroutine htr_lsReInitBlk
