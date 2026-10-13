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
subroutine ib_velGfm3d_fixed(lmda, velx, vely, velz, px, py, pz, &
                             dt, coeff, buffer, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)
   implicit none
   real, dimension(:, :, :), intent(inout) :: velx, vely, velz
   real, dimension(:, :, :), intent(in) :: lmda
   real, dimension(:, :, :), intent(in) :: px, py, pz
   real, intent(in) :: dt, dx, dy, dz, coeff, buffer(3)
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2

   integer :: i, j, k
   real :: lmdax, lmday, lmdaz, weight

   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2+1
            lmdax = (1./2)*(lmda(i, j, k)+lmda(i-1, j, k))
            weight = 1-tanh(lmdax/buffer(1))**2
            velx(i, j, k) = (1-weight)*velx(i, j, k)
         end do
      end do
   end do

   do k = kz1, kz2
      do j = jy1, jy2+1
         do i = ix1, ix2
            lmday = (1./2)*(lmda(i, j, k)+lmda(i, j-1, k))
            weight = 1-tanh(lmday/buffer(2))**2
            vely(i, j, k) = (1-weight)*vely(i, j, k)
         end do
      end do
   end do

   do k = kz1, kz2+1
      do j = jy1, jy2
         do i = ix1, ix2
            lmdaz = (1./2)*(lmda(i, j, k)+lmda(i, j, k-1))
            weight = 1-tanh(lmdaz/buffer(3))**2
            velz(i, j, k) = (1-weight)*velz(i, j, k)
         end do
      end do
   end do

   do k = kz1-1, kz2
      do j = jy1-1, jy2
         do i = ix1-1, ix2

            if (lmda(i, j, k) .lt. 0. .and. lmda(i+1, j, k) .ge. 0.) then
               velx(i+1, j, k) = velx(i+1, j, k)+dt*px(i+1, j, k)
            end if

            if (lmda(i, j, k) .ge. 0. .and. lmda(i+1, j, k) .lt. 0.) then
               velx(i+1, j, k) = velx(i+1, j, k)+dt*px(i+1, j, k)
            end if

            if (lmda(i, j, k) .lt. 0. .and. lmda(i, j+1, k) .ge. 0.) then
               vely(i, j+1, k) = vely(i, j+1, k)+dt*py(i, j+1, k)
            end if

            if (lmda(i, j, k) .ge. 0. .and. lmda(i, j+1, k) .lt. 0.) then
               vely(i, j+1, k) = vely(i, j+1, k)+dt*py(i, j+1, k)
            end if

            if (lmda(i, j, k) .lt. 0. .and. lmda(i, j, k+1) .ge. 0.) then
               velz(i, j, k+1) = velz(i, j, k+1)+dt*pz(i, j, k+1)
            end if

            if (lmda(i, j, k) .ge. 0. .and. lmda(i, j, k+1) .lt. 0.) then
               velz(i, j, k+1) = velz(i, j, k+1)+dt*pz(i, j, k+1)
            end if

         end do
      end do
   end do

end subroutine ib_velGfm3d_fixed
