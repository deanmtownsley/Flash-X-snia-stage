!!****if* source/physics/Multiphase/MultiphaseEvap/mph_evapVelocity
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
!!
!!
!!******
subroutine mph_evapVelocity2d(uni, vni, rhox, rhoy, normx, normy, mflux, ix1, ix2, jy1, jy2)

   !--Arugment List------------------------------
   implicit none
   real, dimension(:, :, :), intent(inout) :: uni, vni
   real, dimension(:, :, :), intent(in)    :: rhox, rhoy
   real, dimension(:, :, :), intent(in)    :: mflux, normx, normy
   integer, intent(in)                   :: ix1, ix2, jy1, jy2

   integer, parameter :: kz1 = 1

   !------U-COMPONENT--------
   uni(ix1:ix2 + 1, jy1:jy2, kz1) = uni(ix1:ix2 + 1, jy1:jy2, kz1) + &
                                    (mflux(ix1 - 1:ix2, jy1:jy2, kz1) + mflux(ix1:ix2 + 1, jy1:jy2, kz1))/2.0d0* &
                                    (normx(ix1 - 1:ix2, jy1:jy2, kz1) + normx(ix1:ix2 + 1, jy1:jy2, kz1))/2.0d0* &
                                    (rhox(ix1:ix2 + 1, jy1:jy2, kz1))

   !------V-COMPONENT--------
   vni(ix1:ix2, jy1:jy2 + 1, kz1) = vni(ix1:ix2, jy1:jy2 + 1, kz1) + &
                                    (mflux(ix1:ix2, jy1 - 1:jy2, kz1) + mflux(ix1:ix2, jy1:jy2 + 1, kz1))/2.0d0* &
                                    (normy(ix1:ix2, jy1 - 1:jy2, kz1) + normy(ix1:ix2, jy1:jy2 + 1, kz1))/2.0d0* &
                                    (rhoy(ix1:ix2, jy1:jy2 + 1, kz1))

end subroutine mph_evapVelocity2d

subroutine mph_evapVelocity3d(uni, vni, wni, rhox, rhoy, rhoz, normx, normy, normz, mflux, ix1, ix2, jy1, jy2, kz1, kz2)

   !--Arugment List------------------------------
   implicit none
   real, dimension(:, :, :), intent(inout) :: uni, vni, wni
   real, dimension(:, :, :), intent(in)    :: rhox, rhoy, rhoz
   real, dimension(:, :, :), intent(in)    :: mflux, normx, normy, normz
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2

   !------U-COMPONENT--------
   uni(ix1:ix2 + 1, jy1:jy2, kz1:kz2) = uni(ix1:ix2 + 1, jy1:jy2, kz1:kz2) + &
                                        (mflux(ix1 - 1:ix2, jy1:jy2, kz1:kz2) + mflux(ix1:ix2 + 1, jy1:jy2, kz1:kz2))/2.0d0* &
                                        (normx(ix1 - 1:ix2, jy1:jy2, kz1:kz2) + normx(ix1:ix2 + 1, jy1:jy2, kz1:kz2))/2.0d0* &
                                        (rhox(ix1:ix2 + 1, jy1:jy2, kz1:kz2))

   !------V-COMPONENT--------
   vni(ix1:ix2, jy1:jy2 + 1, kz1:kz2) = vni(ix1:ix2, jy1:jy2 + 1, kz1:kz2) + &
                                        (mflux(ix1:ix2, jy1 - 1:jy2, kz1:kz2) + mflux(ix1:ix2, jy1:jy2 + 1, kz1:kz2))/2.0d0* &
                                        (normy(ix1:ix2, jy1 - 1:jy2, kz1:kz2) + normy(ix1:ix2, jy1:jy2 + 1, kz1:kz2))/2.0d0* &
                                        (rhoy(ix1:ix2, jy1:jy2 + 1, kz1:kz2))

   !------W-COMPONENT--------
   wni(ix1:ix2, jy1:jy2, kz1:kz2 + 1) = wni(ix1:ix2, jy1:jy2, kz1:kz2 + 1) + &
                                        (mflux(ix1:ix2, jy1:jy2, kz1 - 1:kz2) + mflux(ix1:ix2, jy1:jy2, kz1:kz2 + 1))/2.0d0* &
                                        (normz(ix1:ix2, jy1:jy2, kz1 - 1:kz2) + normz(ix1:ix2, jy1:jy2, kz1:kz2 + 1))/2.0d0* &
                                        (rhoz(ix1:ix2, jy1:jy2, kz1:kz2 + 1))

end subroutine mph_evapVelocity3d
