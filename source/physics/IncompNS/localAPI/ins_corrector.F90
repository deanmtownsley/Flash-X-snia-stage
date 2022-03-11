!!
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
!!
!!**
subroutine ins_corrector_constdens(uni, vni, wni, pxn1, pyn1, pzn1, p, ix1, ix2, jy1, jy2, kz1, kz2, &
                                   dt, dx, dy, dz)
   ! This routine computes the corrected divergence-free velocities.
   implicit none
   INTEGER, INTENT(IN) :: ix1, ix2, jy1, jy2, kz1, kz2
   REAL, INTENT(IN) :: dt, dx, dy, dz
   REAL, DIMENSION(:, :, :), INTENT(IN) :: p
   REAL, DIMENSION(:, :, :), INTENT(INOUT) :: uni, vni, wni, pxn1, pyn1, pzn1
end subroutine ins_corrector_constdens

subroutine ins_corrector_vardens(uni, vni, wni, sigx, sigy, sigz, pxn1, pyn1, pzn1, &
                                 pxn2, pyn2, pzn2, &
                                 rhox, rhoy, rhoz, p, rhoGas, dt, dx, dy, dz, &
                                 ix1, ix2, jy1, jy2, kz1, kz2)

   implicit none
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
   real, intent(in) :: dt, dx, dy, dz, rhoGas
   real, dimension(:, :, :), intent(in) :: p
   real, dimension(:, :, :), intent(in) :: rhox, rhoy, rhoz
   real, dimension(:, :, :), intent(in) :: sigx, sigy, sigz
   real, dimension(:, :, :), intent(inout) :: pxn1, pxn2, pyn1, pyn2, pzn1, pzn2
   real, dimension(:, :, :), intent(inout) :: uni, vni, wni
end subroutine ins_corrector_vardens
