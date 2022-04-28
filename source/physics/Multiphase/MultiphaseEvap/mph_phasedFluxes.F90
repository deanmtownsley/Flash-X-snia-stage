subroutine mph_phasedFluxes(Tn, rhs, phi, dt, ix1, ix2, jy1, jy2, kz1, kz2)
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

   implicit none
   real, dimension(:, :, :), intent(inout) :: Tn
   real, dimension(:, :, :), intent(in) :: rhs, phi
   real, intent(in) :: dt
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2

   real, dimension(ix1:ix2, jy1:jy2, kz1:kz2) :: pf

   pf(ix1:ix2, jy1:jy2, kz1:kz2) = (sign(1.0, phi(ix1:ix2, jy1:jy2, kz1:kz2)) + 1.0)/2.0

   Tn(ix1:ix2, jy1:jy2, kz1:kz2) = Tn(ix1:ix2, jy1:jy2, kz1:kz2) &
                                   + dt*rhs(ix1:ix2, jy1:jy2, kz1:kz2)*pf(ix1:ix2, jy1:jy2, kz1:kz2)

end subroutine mph_phasedFluxes
