!!****if* source/physics/Multiphase/localAPI/mph_evapInterface
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
!! NAME
!!
!!
!! SYNOPSIS
!!  mph_evapInterface()
!!
!! DESCRIPTION
!!  This is an interface specific for the Multiphase
!!  module that defines its private interfaces.
!!
!!***
Module mph_evapInterface

   interface
      subroutine mph_evapInit()
         implicit none
      end subroutine mph_evapInit
   end interface

   interface
      subroutine mph_evapFinalize()
         implicit none
      end subroutine mph_evapFinalize
   end interface

   interface
      subroutine mph_tempGfm2d(phi, nx, ny, Tcoeff, T, Tfrc, Tnl, Tng, Tsat, dx, dy, ix1, ix2, jy1, jy2)
         implicit none
         real, dimension(:, :, :), intent(inout) :: Tfrc, Tnl, Tng
         real, dimension(:, :, :), intent(in) :: phi, T, nx, ny, Tcoeff
         real, intent(in) :: Tsat, dx, dy
         integer, intent(in) :: ix1, ix2, jy1, jy2
      end subroutine mph_tempGfm2d
   end interface

   interface
      subroutine mph_tempGfm3d(phi, nx, ny, nz, Tcoeff, T, Tfrc, Tnl, Tng, Tsat, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)
         implicit none
         real, dimension(:, :, :), intent(inout) :: Tfrc, Tnl, Tng
         real, dimension(:, :, :), intent(in) :: phi, T, nx, ny, nz, Tcoeff
         real, intent(in) :: Tsat, dx, dy, dz
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
      end subroutine mph_tempGfm3d
   end interface

   interface
      subroutine mph_phasedFluxes(Tn, rhs, phi, dt, ix1, ix2, jy1, jy2, kz1, kz2)
         implicit none
         real, dimension(:, :, :), intent(inout) :: Tn
         real, dimension(:, :, :), intent(in) :: rhs, phi
         real, intent(in) :: dt
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
      end subroutine mph_phasedFluxes
   end interface

   interface mph_evapVelocity
      subroutine mph_evapVelocity2d(uni, vni, rhox, rhoy, normx, normy, mflux, ix1, ix2, jy1, jy2)
         implicit none
         real, dimension(:, :, :), intent(inout) :: uni, vni
         real, dimension(:, :, :), intent(in)    :: rhox, rhoy
         real, dimension(:, :, :), intent(in)    :: mflux, normx, normy
         integer, intent(in)                   :: ix1, ix2, jy1, jy2
      end subroutine mph_evapVelocity2d

      subroutine mph_evapVelocity3d(uni, vni, wni, rhox, rhoy, rhoz, normx, normy, normz, mflux, ix1, ix2, jy1, jy2, kz1, kz2)
         implicit none
         real, dimension(:, :, :), intent(inout) :: uni, vni, wni
         real, dimension(:, :, :), intent(in)    :: rhox, rhoy, rhoz
         real, dimension(:, :, :), intent(in)    :: mflux, normx, normy, normz
         integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
      end subroutine mph_evapVelocity3d
   end interface

   interface
      subroutine mph_evapVelForcing2d(uni, vni, rhox, rhoy, visc, normx, normy, mflux, &
                                      ru1, dt, dx, dy, ix1, ix2, jy1, jy2)
         implicit none
         real, dimension(:, :, :), intent(inout) :: uni, vni
         real, dimension(:, :, :), intent(in)    :: rhox, rhoy
         real, dimension(:, :, :), intent(in)    :: visc, normx, normy, mflux
         real                                  :: ru1, dt, dx, dy
         integer, intent(in)                   :: ix1, ix2, jy1, jy2
      end subroutine mph_evapVelForcing2d
   end interface

   interface
      subroutine mph_evapVelForcing3d(uni, vni, wni, rhox, rhoy, rhoz, visc, normx, normy, normz, mflux, &
                                      ru1, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)
         implicit none
         real, dimension(:, :, :), intent(inout) :: uni, vni, wni
         real, dimension(:, :, :), intent(in)    :: rhox, rhoy, rhoz
         real, dimension(:, :, :), intent(in)    :: visc, normx, normy, normz, mflux
         real                                  :: ru1, dt, dx, dy, dz
         integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
      end subroutine mph_evapVelForcing3d
   end interface

   interface mph_evapDivergence
      subroutine mph_evapDivergence2d(divv, rhox, rhoy, normx, normy, mflux, dx, dy, ix1, ix2, jy1, jy2)
         implicit none
         real, dimension(:, :, :), intent(inout) :: divv
         real, dimension(:, :, :), intent(in)    :: rhox, rhoy
         real, dimension(:, :, :), intent(in)    :: mflux, normx, normy
         real, intent(in)                      :: dx, dy
         integer, intent(in)                   :: ix1, ix2, jy1, jy2
      end subroutine mph_evapDivergence2d

      subroutine mph_evapDivergence3d(divv, rhox, rhoy, rhoz, normx, normy, normz, mflux, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)
         implicit none
         real, dimension(:, :, :), intent(inout) :: divv
         real, dimension(:, :, :), intent(in)    :: rhox, rhoy, rhoz
         real, dimension(:, :, :), intent(in)    :: mflux, normx, normy, normz
         real, intent(in)                      :: dx, dy, dz
         integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
      end subroutine mph_evapDivergence3d
   end interface

   interface
      subroutine mph_setEvapJumps2d(phi, sigx, sigy, mflux, rhoGas, dx, dy, ix1, ix2, jy1, jy2)
         implicit none
         integer, intent(in) :: ix1, ix2, jy1, jy2
         real, intent(in) :: dx, dy, rhoGas
         real, dimension(:, :, :), intent(in) :: phi, mflux
         real, dimension(:, :, :), intent(inout) :: sigx, sigy
      end subroutine mph_setEvapJumps2d
   end interface

   interface
      subroutine mph_setEvapJumps3d(phi, sigx, sigy, sigz, mflux, rhoGas, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)
         implicit none
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
         real, intent(in) :: dx, dy, dz, rhoGas
         real, dimension(:, :, :), intent(in) :: phi, mflux
         real, dimension(:, :, :), intent(inout) :: sigx, sigy, sigz
      end subroutine mph_setEvapJumps3d
   end interface

End module mph_evapInterface
