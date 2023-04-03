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
subroutine Stencils_lsCurvature3d(crv, phi, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)
   !
   implicit none

   !-----Argument list-------------------
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
   real, intent(in) :: dx, dy, dz
   real, dimension(:, :, :), intent(in) :: phi
   real, dimension(:, :, :), intent(inout) :: crv

   !-------Local variables---------------
   integer :: i, j, k
   real :: rPhiXN, rPhiXE, rPhiXS, rPhiXW, &
           rPhiYN, rPhiYE, rPhiYS, rPhiYW, &
           rMagN, rMagE, rMagS, rMagW, &
           rPhiZN, rPhiZE, rPhiZS, rPhiZW, &
           rMagF, rMagB, &
           rPhiXF, rPhiXB, rPhiYF, rPhiYB, rPhiZF, rPhiZB
   real, parameter :: eps = 1E-13

   do k = kz1 + 1, kz2 - 1
      do j = jy1 + 1, jy2 - 1
         do i = ix1 + 1, ix2 - 1
            !-------------------------------------------------
            !--3 phi gradients per face method
            !-------------------------------------------------

            !---Compute [d(phi)/dx] on all faces
            rPhiXE = 1./dx*(phi(i + 1, j, k) - phi(i, j, k))
            rPhiXW = 1./dx*(phi(i, j, k) - phi(i - 1, j, k))
            rPhiXN = 1./4./dx*((phi(i + 1, j + 1, k) - phi(i - 1, j + 1, k)) &
                               + (phi(i + 1, j, k) - phi(i - 1, j, k)))
            rPhiXS = 1./4./dx*((phi(i + 1, j, k) - phi(i - 1, j, k)) &
                               + (phi(i + 1, j - 1, k) - phi(i - 1, j - 1, k)))
            rPhiXF = 1./4./dx*((phi(i + 1, j, k + 1) - phi(i - 1, j, k + 1)) &
                               + (phi(i + 1, j, k) - phi(i - 1, j, k)))
            rPhiXB = 1./4./dx*((phi(i + 1, j, k) - phi(i - 1, j, k)) &
                               + (phi(i + 1, j, k - 1) - phi(i - 1, j, k - 1)))

            !---Compute [d(phi)/dy] on all faces
            rPhiYE = 1./4./dy*((phi(i + 1, j + 1, k) - phi(i + 1, j - 1, k)) &
                               + (phi(i, j + 1, k) - phi(i, j - 1, k)))
            rPhiYW = 1./4./dy*((phi(i, j + 1, k) - phi(i, j - 1, k)) &
                               + (phi(i - 1, j + 1, k) - phi(i - 1, j - 1, k)))
            rPhiYN = 1./dy*(phi(i, j + 1, k) - phi(i, j, k))
            rPhiYS = 1./dy*(phi(i, j, k) - phi(i, j - 1, k))
            rPhiYF = 1./4./dy*((phi(i, j + 1, k + 1) - phi(i, j - 1, k + 1)) &
                               + (phi(i, j + 1, k) - phi(i, j - 1, k)))
            rPhiYB = 1./4./dy*((phi(i, j + 1, k) - phi(i, j - 1, k)) &
                               + (phi(i, j + 1, k - 1) - phi(i, j - 1, k - 1)))

            !----Compute [d(phi)/dz] on all faces
            rPhiZE = 1./4./dz*((phi(i + 1, j, k + 1) - phi(i + 1, j, k - 1)) &
                               + (phi(i, j, k + 1) - phi(i, j, k - 1)))
            rPhiZW = 1./4./dz*((phi(i, j, k + 1) - phi(i, j, k - 1)) &
                               + (phi(i - 1, j, k + 1) - phi(i - 1, j, k - 1)))
            rPhiZN = 1./4./dz*((phi(i, j + 1, k + 1) - phi(i, j + 1, k - 1)) &
                               + (phi(i, j, k + 1) - phi(i, j, k - 1)))
            rPhiZS = 1./4./dz*((phi(i, j, k + 1) - phi(i, j, k - 1)) &
                               + (phi(i, j - 1, k + 1) - phi(i, j - 1, k - 1)))
            rPhiZF = 1./dz*(phi(i, j, k + 1) - phi(i, j, k))
            rPhiZB = 1./dz*(phi(i, j, k) - phi(i, j, k - 1))

            !----Compute the magnitude of the normal for ALL faces
            rMagE = sqrt(rPhiXE**2.+rPhiYE**2.+rPhiZE**2.) + eps
            rMagW = sqrt(rPhiXW**2.+rPhiYW**2.+rPhiZW**2.) + eps
            rMagN = sqrt(rPhiXN**2.+rPhiYN**2.+rPhiZN**2.) + eps
            rMagS = sqrt(rPhiXS**2.+rPhiYS**2.+rPhiZS**2.) + eps
            rMagF = sqrt(rPhiXF**2.+rPhiYF**2.+rPhiZF**2.) + eps
            rMagB = sqrt(rPhiXB**2.+rPhiYB**2.+rPhiZB**2.) + eps

            !------------------------------------------------------------
            !--kpd--Finally, compue the curvature, K=grad(s)/||grad(s)||
            crv(i, j, k) = 1./dx*(rPhiXE/rMagE - rPhiXW/rMagW) &
                           + 1./dy*(rPhiYN/rMagN - rPhiYS/rMagS) &
                           + 1./dz*(rPhiZF/rMagF - rPhiZB/rMagB)
            !------------------------------------------------------------
         end do
      end do
   end do

end subroutine Stencils_lsCurvature3d
