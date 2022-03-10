subroutine mph_tempGfm3d(phi, nx, ny, nz, Tcoeff, T, Tfrc, Tnl, Tng, Tsat, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)
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

   implicit none
   real, dimension(:, :, :), intent(inout) :: Tfrc, Tnl, Tng
   real, dimension(:, :, :), intent(in) :: phi, T, nx, ny, nz, Tcoeff
   real, intent(in) :: Tsat, dx, dy, dz
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2

   integer :: i, j, k
   real, parameter :: tol = 0.01
   real :: thxp, thxm, thyp, thym, thzp, thzm
   real :: Txplus, Txmins, Typlus, Tymins, Tzplus, Tzmins
   real :: Tx, Ty, Tz, pf

   do k = kz1 + 1, kz2 - 1
      do j = jy1 + 1, jy2 - 1
         do i = ix1 + 1, ix2 - 1
            thxp = abs(phi(i, j, k))/(abs(phi(i, j, k)) + abs(phi(i + 1, j, k)))
            thxm = abs(phi(i, j, k))/(abs(phi(i, j, k)) + abs(phi(i - 1, j, k)))
            thyp = abs(phi(i, j, k))/(abs(phi(i, j, k)) + abs(phi(i, j + 1, k)))
            thym = abs(phi(i, j, k))/(abs(phi(i, j, k)) + abs(phi(i, j - 1, k)))
            thzp = abs(phi(i, j, k))/(abs(phi(i, j, k)) + abs(phi(i, j, k + 1)))
            thzm = abs(phi(i, j, k))/(abs(phi(i, j, k)) + abs(phi(i, j, k - 1)))

            Txplus = T(i + 1, j, k)
            Txmins = T(i - 1, j, k)
            Typlus = T(i, j + 1, k)
            Tymins = T(i, j - 1, k)
            Tzplus = T(i, j, k + 1)
            Tzmins = T(i, j, k - 1)

            ! Case 1 !
            if (phi(i, j, k)*phi(i + 1, j, k) .le. 0.d0) then
               Txplus = T(i, j, k) + (Tsat - T(i, j, k))/max(tol, thxp)
               Tfrc(i, j, k) = Tfrc(i, j, k) + Tcoeff(i, j, k)*(Txplus - T(i + 1, j, k))/(dx**2)
            end if

            ! Case 2 !
            if (phi(i, j, k)*phi(i - 1, j, k) .le. 0.d0) then
               Txmins = T(i, j, k) + (Tsat - T(i, j, k))/max(tol, thxm)
               Tfrc(i, j, k) = Tfrc(i, j, k) + Tcoeff(i, j, k)*(Txmins - T(i - 1, j, k))/(dx**2)
            end if

            ! Case 3 !
            if (phi(i, j, k)*phi(i, j + 1, k) .le. 0.d0) then
               Typlus = T(i, j, k) + (Tsat - T(i, j, k))/max(tol, thyp)
               Tfrc(i, j, k) = Tfrc(i, j, k) + Tcoeff(i, j, k)*(Typlus - T(i, j + 1, k))/(dy**2)
            end if

            ! Case 4 !
            if (phi(i, j, k)*phi(i, j - 1, k) .le. 0.d0) then
               Tymins = T(i, j, k) + (Tsat - T(i, j, k))/max(tol, thym)
               Tfrc(i, j, k) = Tfrc(i, j, k) + Tcoeff(i, j, k)*(Tymins - T(i, j - 1, k))/(dy**2)
            end if

            ! Case 5 !
            if (phi(i, j, k)*phi(i, j, k + 1) .le. 0.d0) then
               Tzplus = T(i, j, k) + (Tsat - T(i, j, k))/max(tol, thzp)
               Tfrc(i, j, k) = Tfrc(i, j, k) + Tcoeff(i, j, k)*(Tzplus - T(i, j, k + 1))/(dz**2)
            end if

            ! Case 6 !
            if (phi(i, j, k)*phi(i, j, k - 1) .le. 0.d0) then
               Tzmins = T(i, j, k) + (Tsat - T(i, j, k))/max(tol, thzm)
               Tfrc(i, j, k) = Tfrc(i, j, k) + Tcoeff(i, j, k)*(Tzmins - T(i, j, k - 1))/(dz**2)
            end if

            Tx = (Txplus - Txmins)/(2*dx)
            Ty = (Typlus - Tymins)/(2*dy)
            Tz = (Tzplus - Tzmins)/(2*dz)

            pf = (sign(1.0, phi(i, j, k)) + 1.0)/2.0

            Tnl(i, j, k) = (1 - pf)*(nx(i, j, k)*Tx + ny(i, j, k)*Ty + nz(i, j, k)*Tz)
            Tng(i, j, k) = pf*(-nx(i, j, k)*Tx - ny(i, j, k)*Ty - nz(i, j, k)*Tz)

         end do
      end do
   end do

end subroutine mph_tempGfm3d
