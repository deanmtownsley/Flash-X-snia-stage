!!****if* source/physics/Multiphase/MultiphaseEvap/mph_evapVelForcing2d
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
subroutine mph_evapVelForcing2d(uni, vni, rhox, rhoy, visc, normx, normy, mflux, &
                                ru1, dt, dx, dy, ix1, ix2, jy1, jy2)

   implicit none
   real, dimension(:, :, :), intent(inout) :: uni, vni
   real, dimension(:, :, :), intent(in)    :: rhox, rhoy
   real, dimension(:, :, :), intent(in)    :: visc, normx, normy, mflux
   real                                  :: ru1, dt, dx, dy
   integer, intent(in)                   :: ix1, ix2, jy1, jy2

   real :: aicc, aixr, aixl, aiyr, aiyl
   real :: dsdxp, dsdxm, dsdyp, dsdym
   real :: txxp, txxm, tyyp, tyym
   real :: dx1, dy1, Mdens
   integer :: i, j
   integer, parameter :: kz1 = 1

   dx1 = 1./dx
   dy1 = 1./dy

   !------U-COMPONENT--------
   do j = jy1, jy2
      do i = ix1, ix2 + 1

         aicc = rhox(i, j, kz1)* &
                0.5*(mflux(i, j, kz1) + mflux(i - 1, j, kz1))* &
                0.5*(normx(i, j, kz1) + normx(i - 1, j, kz1))

         aixr = rhox(i + 1, j, kz1)* &
                0.5*(mflux(i, j, kz1) + mflux(i - 1, j, kz1))* &
                0.5*(normx(i, j, kz1) + normx(i - 1, j, kz1))

         aixl = rhox(i - 1, j, kz1)* &
                0.5*(mflux(i, j, kz1) + mflux(i - 1, j, kz1))* &
                0.5*(normx(i, j, kz1) + normx(i - 1, j, kz1))

         aiyr = rhox(i, j + 1, kz1)* &
                0.5*(mflux(i, j, kz1) + mflux(i - 1, j, kz1))* &
                0.5*(normx(i, j, kz1) + normx(i - 1, j, kz1))

         aiyl = rhox(i, j - 1, kz1)* &
                0.5*(mflux(i, j, kz1) + mflux(i - 1, j, kz1))* &
                0.5*(normx(i, j, kz1) + normx(i - 1, j, kz1))

         dsdxp = (aixr - aicc)*dx1
         dsdxm = (aicc - aixl)*dx1
         dsdyp = (aiyr - aicc)*dy1
         dsdym = (aicc - aiyl)*dy1

         !-Variable Viscosity Implementation (ru1 is 1/Re)
         txxp = ru1*visc(i, j, kz1)*dsdxp
         txxm = ru1*visc(i - 1, j, kz1)*dsdxm
         tyyp = ru1*0.25*(visc(i, j, kz1) + visc(i - 1, j, kz1) + visc(i - 1, j + 1, kz1) + visc(i, j + 1, kz1))*dsdyp
         tyym = ru1*0.25*(visc(i, j, kz1) + visc(i - 1, j, kz1) + visc(i - 1, j - 1, kz1) + visc(i, j - 1, kz1))*dsdym

         Mdens = rhox(i, j, kz1)

         ! calculate forcing on u-momentum
         uni(i, j, kz1) = uni(i, j, kz1) &
                          + dt*Mdens*(txxp - txxm)*dx1 &
                          + dt*Mdens*(tyyp - tyym)*dy1
      end do
   end do

   !++++++++++  V-COMPONENT  ++++++++++
   do j = jy1, jy2 + 1
      do i = ix1, ix2

         aicc = rhoy(i, j, kz1)* &
                0.5*(mflux(i, j, kz1) + mflux(i, j - 1, kz1))* &
                0.5*(normy(i, j, kz1) + normy(i, j - 1, kz1))

         aixr = rhoy(i + 1, j, kz1)* &
                0.5*(mflux(i, j, kz1) + mflux(i, j - 1, kz1))* &
                0.5*(normy(i, j, kz1) + normy(i, j - 1, kz1))

         aixl = rhoy(i - 1, j, kz1)* &
                0.5*(mflux(i, j, kz1) + mflux(i, j - 1, kz1))* &
                0.5*(normy(i, j, kz1) + normy(i, j - 1, kz1))

         aiyr = rhoy(i, j + 1, kz1)* &
                0.5*(mflux(i, j, kz1) + mflux(i, j - 1, kz1))* &
                0.5*(normy(i, j, kz1) + normy(i, j - 1, kz1))

         aiyl = rhoy(i, j - 1, kz1)* &
                0.5*(mflux(i, j, kz1) + mflux(i, j - 1, kz1))* &
                0.5*(normy(i, j, kz1) + normy(i, j - 1, kz1))

         dsdxp = (aixr - aicc)*dx1
         dsdxm = (aicc - aixl)*dx1
         dsdyp = (aiyr - aicc)*dy1
         dsdym = (aicc - aiyl)*dy1

         !- Variable Viscosity Implementation (ru1 is 1/Re)
         txxp = ru1*0.25*(visc(i, j, kz1) + visc(i + 1, j, kz1) + visc(i, j - 1, kz1) + visc(i + 1, j - 1, kz1))*dsdxp
         txxm = ru1*0.25*(visc(i, j, kz1) + visc(i - 1, j, kz1) + visc(i, j - 1, kz1) + visc(i - 1, j - 1, kz1))*dsdxm
         tyyp = ru1*visc(i, j, kz1)*dsdyp
         tyym = ru1*visc(i, j - 1, kz1)*dsdym

         Mdens = rhoy(i, j, kz1)

         ! calculate forcing on u-momentum
         vni(i, j, kz1) = vni(i, j, kz1) &
                          + dt*Mdens*(txxp - txxm)*dx1 &
                          + dt*Mdens*(tyyp - tyym)*dy1
      end do
   end do

end subroutine mph_evapVelForcing2d
