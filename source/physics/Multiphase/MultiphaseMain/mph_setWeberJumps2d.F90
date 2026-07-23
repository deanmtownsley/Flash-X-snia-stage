!! NOTICE
!!  Copyright 2024 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
subroutine mph_setWeberJumps2d(phi, crv, pf, sigx, sigy, dx, dy, invWbr, rhoGas, ix1, ix2, jy1, jy2, tol)
#include "Simulation.h"
   !
   implicit none

   !-----Argument list-------------------
   integer, intent(in) :: ix1, ix2, jy1, jy2
   real, intent(in) :: dx, dy, invWbr, rhoGas
   real, dimension(:, :, :), intent(in) :: phi, pf
   real, dimension(:, :, :), intent(inout) :: sigx, sigy, crv
   real, intent(in) :: tol

   !-------Local variables---------------
   integer :: icrv(NXB+2*NGUARD,NYB+2*NGUARD,1)
   real :: th, xijl, xijr, &
           cri, xij, yij, yijl, yijr
   integer :: i, j, k
   real, parameter :: eps=1E-13

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
   icrv = 0

   !--Need to loop through one guard cell on each side to set jumps
   !---when they cross block boundaries
   do j = jy1 - 1, jy2
      do i = ix1 - 1, ix2
         !--------------------------------------------------------------
         !- kpd - pf=0 (water) in current cell and pf=1 (air) in cell to right
         !--------------------------------------------------------------
         if (pf(i, j, k) .eq. 0. .and. pf(i + 1, j, k) .eq. 1.) then
            !          = (+)            = (+)           = (-)
            th = max(tol, abs(phi(i + 1, j, k))/(abs(phi(i + 1, j, k)) + abs(phi(i, j, k))))
            xijl = invWbr*crv(i, j, k)                 !- kpd - sigma*K. Used for jump in pressure
            xijr = invWbr*crv(i + 1, j, k)               !- kpd - sigma*K. Used for jump in pressure
            xij = xijl*th + xijr*(1.-th)             !- kpd - Jump in value
            sigx(i + 1, j, k) = sigx(i + 1, j, k) - xij/dx   !- kpd - sigma*K/rho/dx
            icrv(i, j, k) = 1
            icrv(i + 1, j, k) = 1
         end if

         !--------------------------------------------------------------
         !- kpd - pf=1 in current cell and pf=0 in cell to right
         !--------------------------------------------------------------
         if (pf(i, j, k) .eq. 1. .and. pf(i + 1, j, k) .eq. 0.) then
            !
            th = max(tol, abs(phi(i, j, k))/(abs(phi(i, j, k)) + abs(phi(i + 1, j, k))))
            xijl = invWbr*crv(i, j, k)
            xijr = invWbr*crv(i + 1, j, k)
            xij = xijl*(1.-th) + xijr*th
            sigx(i + 1, j, k) = sigx(i + 1, j, k) + xij/dx
            icrv(i, j, k) = 1
            icrv(i + 1, j, k) = 1
         end if

         !--------------------------------------------------------------
         !- kpd - pf=0 in current cell and pf=1 in cell above
         !--------------------------------------------------------------
         if (pf(i, j, k) .eq. 0. .and. pf(i, j + 1, k) .eq. 1.) then
            !
            th = max(tol, abs(phi(i, j + 1, k))/(abs(phi(i, j + 1, k)) + abs(phi(i, j, k))))
            yijl = invWbr*crv(i, j, k)
            yijr = invWbr*crv(i, j + 1, k)
            yij = yijl*th + yijr*(1.-th)
            sigy(i, j + 1, k) = sigy(i, j + 1, k) - yij/dy
            icrv(i, j, k) = 1
            icrv(i, j + 1, k) = 1
         end if

         !--------------------------------------------------------------
         !- kpd - pf=1 in current cell and pf=0 in cell above
         !--------------------------------------------------------------
         if (pf(i, j, k) .eq. 1. .and. pf(i, j + 1, k) .eq. 0.) then
            !
            th = max(tol, abs(phi(i, j, k))/(abs(phi(i, j, k)) + abs(phi(i, j + 1, k))))
            yijl = invWbr*crv(i, j, k)
            yijr = invWbr*crv(i, j + 1, k)
            yij = yijl*(1.-th) + yijr*th
            sigy(i, j + 1, k) = sigy(i, j + 1, k) + yij/dy
            icrv(i, j, k) = 1
            icrv(i, j + 1, k) = 1
         end if
         !--------------------------------------------------------------
         !--------------------------------------------------------------
      end do
   end do

   crv(ix1:ix2, jy1:jy2, 1) = icrv(ix1:ix2, jy1:jy2, 1)*crv(ix1:ix2, jy1:jy2, 1)

end subroutine mph_setWeberJumps2d
