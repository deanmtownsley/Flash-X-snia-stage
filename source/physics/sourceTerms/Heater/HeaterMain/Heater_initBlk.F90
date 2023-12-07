!!***if* source/physics/sourceTerms/Heater/HeaterMain/Heater_initBlk
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
subroutine Heater_initBlk(xcell, ycell, zcell, ix1, ix2, jy1, jy2, kz1, kz2, temp, phi)

   use Heater_type, ONLY: Heater_type_t
   use Heater_data, ONLY: htr_numHeaters, htr_heaterInfo

   implicit none
   real, dimension(:, :, :), intent(inout) :: temp
   real, dimension(:, :, :), intent(inout), optional :: phi
   real, dimension(:), intent(in) :: xcell, ycell, zcell
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2

   type(Heater_type_t), pointer  :: heater
   integer :: i, j, k, htr, isite
   real    :: idfun, iseedY, iseedX, iseedZ, iradius, iheight

   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2
            do htr = 1, htr_numHeaters

               heater => htr_heaterInfo(htr)

               if (present(phi)) then
                  do isite = 1, heater%numSitesAll
                     iheight = heater%radiusInit(isite)*cos(heater%rcdAngle*acos(-1.0)/180)
                     iradius = heater%radiusInit(isite)
                     iseedX = heater%xSiteInit(isite)
                     iseedZ = heater%zSiteInit(isite)
                     iseedY = heater%ySiteInit(isite)+iheight
                     idfun = iradius-sqrt((xcell(i)-iseedX)**2+(ycell(j)-iseedY)**2+(zcell(k)-iseedZ)**2)
                     phi(i, j, k) = max(phi(i, j, k), idfun)
                  end do
               end if

               if (xcell(i) .ge. heater%xMin .and. &
                   xcell(i) .le. heater%xMax .and. &
                   ycell(j) .le. 0.2 .and. &
                   zcell(k) .ge. heater%zMin .and. &
                   zcell(k) .le. heater%zMax) temp(i, j, k) = (0.2-ycell(j))/0.2

            end do
         end do
      end do
   end do

   return

end subroutine Heater_initBlk
