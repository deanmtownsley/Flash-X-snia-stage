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
subroutine mph_tempGfm2d(phi, nx, ny, Tcoeff, T, Tfrc, Tnl, Tng, Tsat, dx, dy, ix1, ix2, jy1, jy2, tol)
   implicit none
   real, dimension(:, :, :), intent(inout) :: Tfrc, Tnl, Tng
   real, dimension(:, :, :), intent(in) :: phi, T, nx, ny, Tcoeff
   real, intent(in) :: Tsat, dx, dy
   integer, intent(in) :: ix1, ix2, jy1, jy2
   real, intent(in) :: tol
end subroutine mph_tempGfm2d
