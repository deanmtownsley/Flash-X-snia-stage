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
end subroutine ib_velGfm3d_fixed
