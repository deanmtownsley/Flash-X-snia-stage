!!****if* source/physics/Multiphase/localAPI/mph_evapVelocity
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
subroutine mph_evapVelocity2d(uni, vni, rhoc, normx, normy, mflux, ix1, ix2, jy1, jy2)
   implicit none
   real, dimension(:, :, :), intent(inout) :: uni, vni
   real, dimension(:, :, :), intent(in)    :: rhoc
   real, dimension(:, :, :), intent(in)    :: mflux, normx, normy
   integer, intent(in)                   :: ix1, ix2, jy1, jy2
end subroutine mph_evapVelocity2d

subroutine mph_evapVelocity3d(uni, vni, wni, rhoc, normx, normy, normz, mflux, ix1, ix2, jy1, jy2, kz1, kz2)
   implicit none
   real, dimension(:, :, :), intent(inout) :: uni, vni, wni
   real, dimension(:, :, :), intent(in)    :: rhoc
   real, dimension(:, :, :), intent(in)    :: mflux, normx, normy, normz
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
end subroutine mph_evapVelocity3d
