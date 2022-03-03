subroutine mph_setWeberJumps2d(phi, sigx, sigy, dx, dy, invWbr, ix1, ix2, jy1, jy2)
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
   real, intent(in) :: dx, dy, invWbr
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:, :, :), intent(inout) :: sigx, sigy

   !-------Local variables---------------
   real, dimension(ix1:ix2, jy1:jy2, 1) :: crv, pf
   real :: th, aa, xijl, xijr, &
           cri, xij, yij, yijl, yijr
   integer :: i, j, k
   real, parameter :: eps = 1E-13
   real :: rPhiXN, rPhiXE, rPhiXS, rPhiXW, &
           rPhiYN, rPhiYE, rPhiYS, rPhiYW, &
           rMagN, rMagE, rMagS, rMagW

   !--------------------------------------------
   !----------------jump conditions ------------
   !--------------------------------------------
   !xij: jump in value
   !xid: jump in gradient
   !l,r, values at pts left and right of the interface
   !crv: curvature
   !cri: curvature at interface

   !left interface between i and i+1
   !  phase 2 at i, phase 1 at i+1
   !  theta = x_i+1 - x_I

   !right interface between i and i+1
   !  phase 1 at i, phase 2 at i+1
   !  theta = x_I - x_i
   !--------------------------------------------
   !--------------------------------------------

   k = 1
   crv = 0.

   pf(ix1:ix2, jy1:jy2, k) = (sign(1.0, phi(ix1:ix2, jy1:jy2, k)) + 1.0)/2.0

   do j = jy1 + 1, jy2 - 1
      do i = ix1 + 1, ix2 - 1
         !----------------------------------------------------
         !--------------2 phi gradients per face method------
         !----------------------------------------------------
         !        X - Location
         rPhiXE = 1./dx*(phi(i + 1, j, k) - phi(i, j, k))
         rPhiXW = 1./dx*(phi(i, j, k) - phi(i - 1, j, k))
         rPhiXN = 1./4./dx*((phi(i + 1, j + 1, k) - phi(i - 1, j + 1, k)) &
                            + (phi(i + 1, j, k) - phi(i - 1, j, k)))
         rPhiXS = 1./4./dx*((phi(i + 1, j, k) - phi(i - 1, j, k)) &
                            + (phi(i + 1, j - 1, k) - phi(i - 1, j - 1, k)))
         !        Y - Location
         rPhiYN = 1./dy*(phi(i, j + 1, k) - phi(i, j, k))
         rPhiYS = 1./dy*(phi(i, j, k) - phi(i, j - 1, k))
         rPhiYE = 1./4./dy*((phi(i + 1, j + 1, k) - phi(i + 1, j - 1, k)) &
                            + (phi(i, j + 1, k) - phi(i, j - 1, k)))
         rPhiYW = 1./4./dy*((phi(i, j + 1, k) - phi(i, j - 1, k)) &
                            + (phi(i - 1, j + 1, k) - phi(i - 1, j - 1, k)))
         !----------------------------------------------------

         !----Compute the magnitude of the gradient at each face
         rMagE = sqrt(rPhiXE**2.+rPhiYE**2.) + eps
         rMagW = sqrt(rPhiXW**2.+rPhiYW**2.) + eps
         rMagN = sqrt(rPhiXN**2.+rPhiYN**2.) + eps
         rMagS = sqrt(rPhiXS**2.+rPhiYS**2.) + eps

         crv(i, j, k) = 1./dx*(rPhiXE/rMagE - rPhiXW/rMagW) &
                        + 1./dy*(rPhiYN/rMagN - rPhiYS/rMagS)
         !----------------------------------------------------
      end do
   end do

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
            xijl = invWbr*crv(i, j, k)                 !- kpd - sigma*K. Used for jump in pressure
            xijr = invWbr*crv(i + 1, j, k)               !- kpd - sigma*K. Used for jump in pressure
            xij = xijl*th + xijr*(1.-th)             !- kpd - Jump in value
            sigx(i + 1, j, k) = sigx(i + 1, j, k) - xij/dx   !- kpd - sigma*K/rho/dx
         end if

         !--------------------------------------------------------------
         !- kpd - pf=1 in current cell and pf=0 in cell to right
         !--------------------------------------------------------------
         if (pf(i, j, k) .eq. 1. .and. pf(i + 1, j, k) .eq. 0.) then
            !
            th = abs(phi(i, j, k))/(abs(phi(i, j, k)) + abs(phi(i + 1, j, k)))
            xijl = invWbr*crv(i, j, k)
            xijr = invWbr*crv(i + 1, j, k)
            xij = xijl*(1.-th) + xijr*th
            sigx(i + 1, j, k) = sigx(i + 1, j, k) + xij/dx
         end if

         !--------------------------------------------------------------
         !- kpd - pf=0 in current cell and pf=1 in cell above
         !--------------------------------------------------------------
         if (pf(i, j, k) .eq. 0. .and. pf(i, j + 1, k) .eq. 1.) then
            !
            th = abs(phi(i, j + 1, k))/(abs(phi(i, j + 1, k)) + abs(phi(i, j, k)))
            yijl = invWbr*crv(i, j, k)
            yijr = invWbr*crv(i, j + 1, k)
            yij = yijl*th + yijr*(1.-th)
            sigy(i, j + 1, k) = sigy(i, j + 1, k) - yij/dy
         end if

         !--------------------------------------------------------------
         !- kpd - pf=1 in current cell and pf=0 in cell above
         !--------------------------------------------------------------
         if (pf(i, j, k) .eq. 1. .and. pf(i, j + 1, k) .eq. 0.) then
            !
            th = abs(phi(i, j, k))/(abs(phi(i, j, k)) + abs(phi(i, j + 1, k)))
            yijl = invWbr*crv(i, j, k)
            yijr = invWbr*crv(i, j + 1, k)
            yij = yijl*(1.-th) + yijr*th
            sigy(i, j + 1, k) = sigy(i, j + 1, k) + yij/dy
         end if
         !--------------------------------------------------------------
         !--------------------------------------------------------------
      end do
   end do
end subroutine mph_setWeberJumps2d
