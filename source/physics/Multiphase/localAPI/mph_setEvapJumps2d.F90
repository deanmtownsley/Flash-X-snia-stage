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
subroutine mph_setEvapJumps2d(phi, pf, sigx, sigy, mflux, rhoGas, dx, dy, ix1, ix2, jy1, jy2, tol)
   implicit none
   integer, intent(in) :: ix1, ix2, jy1, jy2
   real, intent(in) :: dx, dy, rhoGas
   real, dimension(:, :, :), intent(in) :: phi, mflux, pf
   real, dimension(:, :, :), intent(inout) :: sigx, sigy
   real, intent(in) :: tol
end subroutine mph_setEvapJumps2d
