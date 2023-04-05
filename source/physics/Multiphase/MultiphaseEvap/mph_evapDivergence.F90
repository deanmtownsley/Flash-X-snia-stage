!!****if* source/physics/Multiphase/MultiphaseEvap/mph_evapDivergence
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
!!
!!******
subroutine mph_evapDivergence2d(divv, rhoc, normx, normy, mflux, dx, dy, ix1, ix2, jy1, jy2)

   implicit none
   real, dimension(:, :, :), intent(inout) :: divv
   real, dimension(:, :, :), intent(in)    :: rhoc
   real, dimension(:, :, :), intent(in)    :: mflux, normx, normy
   real, intent(in)                      :: dx, dy
   integer, intent(in)                   :: ix1, ix2, jy1, jy2

   integer :: i, j
   integer, parameter :: k = 1
   real :: rhoxr, rhoxl, rhoyr, rhoyl
   real :: aicx, aicy

   do j = jy1, jy2
      do i = ix1, ix2

         rhoxr = (rhoc(i, j, k) + rhoc(i + 1, j, k))/2.0d0
         rhoxl = (rhoc(i, j, k) + rhoc(i - 1, j, k))/2.0d0
         rhoyr = (rhoc(i, j, k) + rhoc(i, j + 1, k))/2.0d0
         rhoyl = (rhoc(i, j, k) + rhoc(i, j - 1, k))/2.0d0

         aicx = mflux(i, j, k)*normx(i, j, k)
         aicy = mflux(i, j, k)*normy(i, j, k)

         divv(i, j, k) = divv(i, j, k) + aicx*(rhoxr - rhoxl)/dx + aicy*(rhoyr - rhoyl)/dy

      end do
   end do

end subroutine mph_evapDivergence2d

subroutine mph_evapDivergence3d(divv, rhoc, normx, normy, normz, mflux, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)

   implicit none
   real, dimension(:, :, :), intent(inout) :: divv
   real, dimension(:, :, :), intent(in)    :: rhoc
   real, dimension(:, :, :), intent(in)    :: mflux, normx, normy, normz
   real, intent(in)                      :: dx, dy, dz
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2

   integer :: i, j
   integer :: k
   real :: rhoxr, rhoxl, rhoyr, rhoyl, rhozr, rhozl
   real :: aicx, aicy, aicz

   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2

            rhoxr = (rhoc(i, j, k) + rhoc(i + 1, j, k))/2.0d0
            rhoxl = (rhoc(i, j, k) + rhoc(i - 1, j, k))/2.0d0
            rhoyr = (rhoc(i, j, k) + rhoc(i, j + 1, k))/2.0d0
            rhoyl = (rhoc(i, j, k) + rhoc(i, j - 1, k))/2.0d0
            rhozr = (rhoc(i, j, k) + rhoc(i, j, k + 1))/2.0d0
            rhozl = (rhoc(i, j, k) + rhoc(i, j, k - 1))/2.0d0

            aicx = mflux(i, j, k)*normx(i, j, k)
            aicy = mflux(i, j, k)*normy(i, j, k)
            aicz = mflux(i, j, k)*normz(i, j, k)

            divv(i, j, k) = divv(i, j, k) + aicx*(rhoxr - rhoxl)/dx + aicy*(rhoyr - rhoyl)/dy + aicz*(rhozr - rhozl)/dz

         end do
      end do
   end do

end subroutine mph_evapDivergence3d
