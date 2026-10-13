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
subroutine Stencils_lsCurvature3d(crv, phi, nrmx, nrmy, nrmz, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)
   !
   implicit none

   !-----Argument list-------------------
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
   real, intent(in) :: dx, dy, dz
   real, dimension(:, :, :), intent(in) :: phi, nrmx, nrmy, nrmz
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

   do k = kz1+1, kz2-1
      do j = jy1+1, jy2-1
         do i = ix1+1, ix2-1

            crv(i, j, k) = (nrmx(i+1, j, k)-nrmx(i-1, j, k))/(2*dx)+ &
                           (nrmy(i, j+1, k)-nrmy(i, j-1, k))/(2*dy)+ &
                           (nrmz(i, j, k+1)-nrmz(i, j, k-1))/(2*dz)

         end do
      end do
   end do

end subroutine Stencils_lsCurvature3d
