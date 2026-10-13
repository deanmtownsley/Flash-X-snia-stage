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
#include "Simulation.h"
#include "constants.h"

subroutine mph_evapVelocity2d(uni, vni, rhoc, normx, normy, mflux, lo, hi)
  
  !--Arugment List------------------------------
  implicit none
  real, dimension(:, :, :), intent(inout) :: uni, vni
  real, dimension(:, :, :), intent(in)    :: rhoc
  real, dimension(:, :, :), intent(in)    :: mflux, normx, normy
  integer, dimension(MDIM), intent(in)    :: lo, hi
  
  integer, parameter :: kz1 = 1
  
  uni(lo(IAXIS):hi(IAXIS) + 1, lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)) = &
       uni(lo(IAXIS):hi(IAXIS) + 1, lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
       (mflux(lo(IAXIS) - 1:hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
       mflux(lo(IAXIS):hi(IAXIS) + 1, lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)))/2.0d0* &
       (normx(lo(IAXIS) - 1:hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
       normx(lo(IAXIS):hi(IAXIS) + 1, lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)))/2.0d0* &
       (rhoc(lo(IAXIS) - 1:hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
       rhoc(lo(IAXIS):hi(IAXIS) + 1, lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)))/2.0d0
  
  !------V-COMPONENT--------
  vni(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS) + 1, lo(KAXIS):hi(KAXIS)) = &
       vni(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS) + 1, lo(KAXIS):hi(KAXIS)) + &
       (mflux(lo(IAXIS):hi(IAXIS), lo(JAXIS) - 1:hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
       mflux(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS) + 1, lo(KAXIS):hi(KAXIS)))/2.0d0* &
       (normy(lo(IAXIS):hi(IAXIS), lo(JAXIS) - 1:hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
       normy(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS) + 1, lo(KAXIS):hi(KAXIS)))/2.0d0* &
       (rhoc(lo(IAXIS):hi(IAXIS), lo(JAXIS) - 1:hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
       rhoc(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS) + 1, lo(KAXIS):hi(KAXIS)))/2.0d0
end subroutine mph_evapVelocity2d

subroutine mph_evapVelocity3d(uni, vni, wni, rhoc, normx, normy, normz, mflux, lo, hi)

   !--Arugment List------------------------------
   implicit none
   real, dimension(:, :, :), intent(inout) :: uni, vni, wni
   real, dimension(:, :, :), intent(in)    :: rhoc
   real, dimension(:, :, :), intent(in)    :: mflux, normx, normy, normz
   integer, dimension(MDIM), intent(in)                   :: lo, hi

   !------U-COMPONENT--------
   uni(lo(IAXIS):hi(IAXIS) + 1, lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)) = &
        uni(lo(IAXIS):hi(IAXIS) + 1, lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
        (mflux(lo(IAXIS) - 1:hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
        mflux(lo(IAXIS):hi(IAXIS) + 1, lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)))/2.0d0* &
        (normx(lo(IAXIS) - 1:hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
        normx(lo(IAXIS):hi(IAXIS) + 1, lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)))/2.0d0* &
        (rhoc(lo(IAXIS) - 1:hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
        rhoc(lo(IAXIS):hi(IAXIS) + 1, lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)))/2.0d0

   !------V-COMPONENT--------
   vni(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS) + 1, lo(KAXIS):hi(KAXIS)) = &
        vni(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS) + 1, lo(KAXIS):hi(KAXIS)) + &
        (mflux(lo(IAXIS):hi(IAXIS), lo(JAXIS) - 1:hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
        mflux(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS) + 1, lo(KAXIS):hi(KAXIS)))/2.0d0* &
        (normy(lo(IAXIS):hi(IAXIS), lo(JAXIS) - 1:hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
        normy(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS) + 1, lo(KAXIS):hi(KAXIS)))/2.0d0* &
        (rhoc(lo(IAXIS):hi(IAXIS), lo(JAXIS) - 1:hi(JAXIS), lo(KAXIS):hi(KAXIS)) + &
        rhoc(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS) + 1, lo(KAXIS):hi(KAXIS)))/2.0d0

   !------W-COMPONENT--------
   wni(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS) + 1) = &
        wni(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS) + 1) + &
        (mflux(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS) - 1:hi(KAXIS)) + &
        mflux(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS) + 1))/2.0d0* &
        (normz(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS) - 1:hi(KAXIS)) + &
        normz(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS) + 1))/2.0d0* &
        (rhoc(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS) - 1:hi(KAXIS)) + &
        rhoc(lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS) + 1))/2.0d0

 end subroutine mph_evapVelocity3d
