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
subroutine Stencils_lsCurvature2d(crv, phi, dx, dy, ix1, ix2, jy1, jy2)
   implicit none

   !-----Argument list-------------------
   integer, intent(in) :: ix1, ix2, jy1, jy2
   real, intent(in) :: dx, dy
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:, :, :), intent(inout) :: crv

   !-------Local variables---------------
   integer :: i, j, k
   real, parameter :: eps = 1E-13
   real :: rPhiXN, rPhiXE, rPhiXS, rPhiXW, &
           rPhiYN, rPhiYE, rPhiYS, rPhiYW, &
           rMagN, rMagE, rMagS, rMagW

   k = 1

   do j = jy1 + 1, jy2 - 1
      do i = ix1 + 1, ix2 - 1
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

end subroutine Stencils_lsCurvature2d
