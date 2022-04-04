subroutine mph_setEvapJumps2d(phi, sigx, sigy, mflux, rhoGas, dx, dy, ix1, ix2, jy1, jy2)
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
   !
   implicit none

   !-----Argument list-------------------
   integer, intent(in) :: ix1, ix2, jy1, jy2
   real, intent(in) :: dx, dy, rhoGas
   real, dimension(:, :, :), intent(in) :: phi, mflux
   real, dimension(:, :, :), intent(inout) :: sigx, sigy

   !-------Local variables---------------
   real, dimension(ix1:ix2, jy1:jy2, 1) :: pf
   real :: th, xijl, xijr, xij, yij, yijl, yijr
   integer :: i, j, k
   real :: bb

   bb = (1./rhoGas) - 1
   k = 1
   pf(ix1:ix2, jy1:jy2, k) = (sign(1.0, phi(ix1:ix2, jy1:jy2, k)) + 1.0)/2.0

   !--Need to loop through one guard cell on each side to set jumps
   !---when they cross block boundaries
   do j = jy1 + 1, jy2 - 2
      do i = ix1 + 1, ix2 - 2
         !--------------------------------------------------------------
         !- kpd - pf=0 (water) in current cell and pf=1 (air) in cell to right
         !--------------------------------------------------------------
         if (pf(i, j, k) .eq. 0. .and. pf(i + 1, j, k) .eq. 1.) then
            !          = (+)            = (+)           = (-)
            th = abs(phi(i + 1, j, k))/(abs(phi(i + 1, j, k)) + abs(phi(i, j, k)))
            xijl = -bb*mflux(i, j, k)*mflux(i, j, k)
            xijr = -bb*mflux(i + 1, j, k)*mflux(i + 1, j, k)
            xij = xijl*th + xijr*(1.-th)
            sigx(i + 1, j, k) = sigx(i + 1, j, k) - xij/dx
         end if

         !--------------------------------------------------------------
         !- kpd - pf=1 in current cell and pf=0 in cell to right
         !--------------------------------------------------------------
         if (pf(i, j, k) .eq. 1. .and. pf(i + 1, j, k) .eq. 0.) then
            !
            th = abs(phi(i, j, k))/(abs(phi(i, j, k)) + abs(phi(i + 1, j, k)))
            xijl = -bb*mflux(i, j, k)*mflux(i, j, k)
            xijr = -bb*mflux(i + 1, j, k)*mflux(i + 1, j, k)
            xij = xijl*(1.-th) + xijr*th
            sigx(i + 1, j, k) = sigx(i + 1, j, k) + xij/dx

         end if

         !--------------------------------------------------------------
         !- kpd - pf=0 in current cell and pf=1 in cell above
         !--------------------------------------------------------------
         if (pf(i, j, k) .eq. 0. .and. pf(i, j + 1, k) .eq. 1.) then
            !
            th = abs(phi(i, j + 1, k))/(abs(phi(i, j + 1, k)) + abs(phi(i, j, k)))
            yijl = -bb*mflux(i, j, k)*mflux(i, j, k)
            yijr = -bb*mflux(i, j + 1, k)*mflux(i, j + 1, k)
            yij = yijl*th + yijr*(1.-th)
            sigy(i, j + 1, k) = sigy(i, j + 1, k) - yij/dy
         end if

         !--------------------------------------------------------------
         !- kpd - pf=1 in current cell and pf=0 in cell above
         !--------------------------------------------------------------
         if (pf(i, j, k) .eq. 1. .and. pf(i, j + 1, k) .eq. 0.) then
            !
            th = abs(phi(i, j, k))/(abs(phi(i, j, k)) + abs(phi(i, j + 1, k)))
            yijl = -bb*mflux(i, j, k)*mflux(i, j, k)
            yijr = -bb*mflux(i, j + 1, k)*mflux(i, j + 1, k)
            yij = yijl*(1.-th) + yijr*th
            sigy(i, j + 1, k) = sigy(i, j + 1, k) + yij/dy
         end if
         !--------------------------------------------------------------
         !--------------------------------------------------------------
      end do
   end do
end subroutine mph_setEvapJumps2d
