!!****if* source/physics/Multiphase/MultiphaseEvap/mph_evapVelForcing3d
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
subroutine mph_evapVelForcing3d(uni, vni, wni, rhox, rhoy, rhoz, visc, normx, normy, normz, mflux, &
                                ru1, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)

   implicit none
   real, dimension(:, :, :), intent(inout) :: uni, vni, wni
   real, dimension(:, :, :), intent(in)    :: rhox, rhoy, rhoz
   real, dimension(:, :, :), intent(in)    :: visc, normx, normy, normz, mflux
   real                                  :: ru1, dt, dx, dy, dz
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2

   real :: aicc, aixr, aixl, aiyr, aiyl, aizr, aizl
   real :: dsdxp, dsdxm, dsdyp, dsdym, dsdzp, dsdzm
   real :: txxp, txxm, tyyp, tyym, tzzp, tzzm
   real :: dx1, dy1, dz1, Mdens
   integer :: i, j, k

   dz1 = 1./dz
   dx1 = 1./dx
   dy1 = 1./dy

   !------U-COMPONENT--------
   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2 + 1

            aicc = rhox(i, j, k)* &
                   0.5*(mflux(i, j, k) + mflux(i - 1, j, k))* &
                   0.5*(normx(i, j, k) + normx(i - 1, j, k))

            aixr = rhox(i + 1, j, k)* &
                   0.5*(mflux(i, j, k) + mflux(i - 1, j, k))* &
                   0.5*(normx(i, j, k) + normx(i - 1, j, k))

            aixl = rhox(i - 1, j, k)* &
                   0.5*(mflux(i, j, k) + mflux(i - 1, j, k))* &
                   0.5*(normx(i, j, k) + normx(i - 1, j, k))

            aiyr = rhox(i, j + 1, k)* &
                   0.5*(mflux(i, j, k) + mflux(i - 1, j, k))* &
                   0.5*(normx(i, j, k) + normx(i - 1, j, k))

            aiyl = rhox(i, j - 1, k)* &
                   0.5*(mflux(i, j, k) + mflux(i - 1, j, k))* &
                   0.5*(normx(i, j, k) + normx(i - 1, j, k))

            aizr = rhox(i, j, k + 1)* &
                   0.5*(mflux(i, j, k) + mflux(i - 1, j, k))* &
                   0.5*(normx(i, j, k) + normx(i - 1, j, k))

            aizl = rhox(i, j, k - 1)* &
                   0.5*(mflux(i, j, k) + mflux(i - 1, j, k))* &
                   0.5*(normx(i, j, k) + normx(i - 1, j, k))

            dsdxp = (aixr - aicc)*dx1
            dsdxm = (aicc - aixl)*dx1
            dsdyp = (aiyr - aicc)*dy1
            dsdym = (aicc - aiyl)*dy1
            dsdzp = (aizr - aicc)*dz1
            dsdzm = (aicc - aizl)*dz1

            !-Variable Viscosity Implementation (ru1 is 1/Re)
            txxp = ru1*visc(i, j, k)*dsdxp
            txxm = ru1*visc(i - 1, j, k)*dsdxm
            tyyp = ru1*0.25*(visc(i, j, k) + visc(i - 1, j, k) + visc(i - 1, j + 1, k) + visc(i, j + 1, k))*dsdyp
            tyym = ru1*0.25*(visc(i, j, k) + visc(i - 1, j, k) + visc(i - 1, j - 1, k) + visc(i, j - 1, k))*dsdym
            tzzp = ru1*0.25*(visc(i, j, k) + visc(i - 1, j, k) + visc(i - 1, j, k + 1) + visc(i, j, k + 1))*dsdzp
            tzzm = ru1*0.25*(visc(i, j, k) + visc(i - 1, j, k) + visc(i - 1, j, k - 1) + visc(i, j, k - 1))*dsdzm

            Mdens = rhox(i, j, k)

            ! calculate forcing on u-momentum
            uni(i, j, k) = uni(i, j, k) &
                           + dt*Mdens*(txxp - txxm)*dx1 &
                           + dt*Mdens*(tyyp - tyym)*dy1 &
                           + dt*Mdens*(tzzp - tzzm)*dz1
         end do
      end do
   end do

   !++++++++++  V-COMPONENT  ++++++++++
   do k = kz1, kz2
      do j = jy1, jy2 + 1
         do i = ix1, ix2

            aicc = rhoy(i, j, k)* &
                   0.5*(mflux(i, j, k) + mflux(i, j - 1, k))* &
                   0.5*(normy(i, j, k) + normy(i, j - 1, k))

            aixr = rhoy(i + 1, j, k)* &
                   0.5*(mflux(i, j, k) + mflux(i, j - 1, k))* &
                   0.5*(normy(i, j, k) + normy(i, j - 1, k))

            aixl = rhoy(i - 1, j, k)* &
                   0.5*(mflux(i, j, k) + mflux(i, j - 1, k))* &
                   0.5*(normy(i, j, k) + normy(i, j - 1, k))

            aiyr = rhoy(i, j + 1, k)* &
                   0.5*(mflux(i, j, k) + mflux(i, j - 1, k))* &
                   0.5*(normy(i, j, k) + normy(i, j - 1, k))

            aiyl = rhoy(i, j - 1, k)* &
                   0.5*(mflux(i, j, k) + mflux(i, j - 1, k))* &
                   0.5*(normy(i, j, k) + normy(i, j - 1, k))

            aizr = rhoy(i, j, k + 1)* &
                   0.5*(mflux(i, j, k) + mflux(i, j - 1, k))* &
                   0.5*(normy(i, j, k) + normy(i, j - 1, k))

            aizl = rhoy(i, j, k - 1)* &
                   0.5*(mflux(i, j, k) + mflux(i, j - 1, k))* &
                   0.5*(normy(i, j, k) + normy(i, j - 1, k))

            dsdxp = (aixr - aicc)*dx1
            dsdxm = (aicc - aixl)*dx1
            dsdyp = (aiyr - aicc)*dy1
            dsdym = (aicc - aiyl)*dy1
            dsdzp = (aizr - aicc)*dz1
            dsdzm = (aicc - aizl)*dz1

            !- Variable Viscosity Implementation (ru1 is 1/Re)
            txxp = ru1*0.25*(visc(i, j, k) + visc(i + 1, j, k) + visc(i, j - 1, k) + visc(i + 1, j - 1, k))*dsdxp
            txxm = ru1*0.25*(visc(i, j, k) + visc(i - 1, j, k) + visc(i, j - 1, k) + visc(i - 1, j - 1, k))*dsdxm
            tyyp = ru1*visc(i, j, k)*dsdyp
            tyym = ru1*visc(i, j - 1, k)*dsdym
            tzzp = ru1*0.25*(visc(i, j, k) + visc(i, j, k + 1) + visc(i, j - 1, k) + visc(i, j - 1, k + 1))*dsdzp
            tzzm = ru1*0.25*(visc(i, j, k) + visc(i, j, k - 1) + visc(i, j - 1, k) + visc(i, j - 1, k - 1))*dsdzm

            Mdens = rhoy(i, j, k)

            ! calculate forcing on u-momentum
            vni(i, j, k) = vni(i, j, k) &
                           + dt*Mdens*(txxp - txxm)*dx1 &
                           + dt*Mdens*(tyyp - tyym)*dy1 &
                           + dt*Mdens*(tzzp - tzzm)*dz1
         end do
      end do
   end do

   !++++++++++  W-COMPONENT  ++++++++++
   do k = kz1, kz2 + 1
      do j = jy1, jy2
         do i = ix1, ix2

            aicc = rhoz(i, j, k)* &
                   0.5*(mflux(i, j, k) + mflux(i, j, k - 1))* &
                   0.5*(normz(i, j, k) + normz(i, j, k - 1))

            aixr = rhoz(i + 1, j, k)* &
                   0.5*(mflux(i, j, k) + mflux(i, j, k - 1))* &
                   0.5*(normz(i, j, k) + normz(i, j, k - 1))

            aixl = rhoz(i - 1, j, k)* &
                   0.5*(mflux(i, j, k) + mflux(i, j, k - 1))* &
                   0.5*(normz(i, j, k) + normz(i, j, k - 1))

            aiyr = rhoz(i, j + 1, k)* &
                   0.5*(mflux(i, j, k) + mflux(i, j, k - 1))* &
                   0.5*(normz(i, j, k) + normz(i, j, k - 1))

            aiyl = rhoz(i, j - 1, k)* &
                   0.5*(mflux(i, j, k) + mflux(i, j, k - 1))* &
                   0.5*(normz(i, j, k) + normz(i, j, k - 1))

            aizr = rhoz(i, j, k + 1)* &
                   0.5*(mflux(i, j, k) + mflux(i, j, k - 1))* &
                   0.5*(normz(i, j, k) + normz(i, j, k - 1))

            aizl = rhoz(i, j, k - 1)* &
                   0.5*(mflux(i, j, k) + mflux(i, j, k - 1))* &
                   0.5*(normz(i, j, k) + normz(i, j, k - 1))

            dsdxp = (aixr - aicc)*dx1
            dsdxm = (aicc - aixl)*dx1
            dsdyp = (aiyr - aicc)*dy1
            dsdym = (aicc - aiyl)*dy1
            dsdzp = (aizr - aicc)*dz1
            dsdzm = (aicc - aizl)*dz1

            !- Variable Viscosity Implementation (ru1 is 1/Re)
            txxp = ru1*0.25*(visc(i, j, k) + visc(i + 1, j, k) + visc(i, j, k - 1) + visc(i + 1, j, k - 1))*dsdxp
            txxm = ru1*0.25*(visc(i, j, k) + visc(i - 1, j, k) + visc(i, j, k - 1) + visc(i - 1, j, k - 1))*dsdxm
            tyyp = ru1*0.25*(visc(i, j, k) + visc(i, j, k - 1) + visc(i, j + 1, k - 1) + visc(i, j + 1, k))*dsdyp
            tyym = ru1*0.25*(visc(i, j, k) + visc(i, j, k - 1) + visc(i, j - 1, k - 1) + visc(i, j - 1, k))*dsdym
            tzzp = ru1*visc(i, j, k)*dsdzp
            tzzm = ru1*visc(i, j, k - 1)*dsdzm

            Mdens = rhoz(i, j, k)

            ! calculate forcing on u-momentum
            wni(i, j, k) = wni(i, j, k) &
                           + dt*Mdens*(txxp - txxm)*dx1 &
                           + dt*Mdens*(tyyp - tyym)*dy1 &
                           + dt*Mdens*(tzzp - tzzm)*dz1
         end do
      end do
   end do

end subroutine mph_evapVelForcing3d
