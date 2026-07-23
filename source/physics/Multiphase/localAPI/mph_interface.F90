!!****if* source/physics/Multiphase/localAPI/mph_interface
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
!!  mph_interface()
!!
!! DESCRIPTION
!!  This is an interface specific for the Multiphase
!!  module that defines its private interfaces.
!!
!!***
Module mph_interface

#include "constants.h"
#include "Simulation.h"
#include "Multiphase.h"

   interface
      subroutine mph_setWeberJumps2d(phi, crv, pf, sigx, sigy, dx, dy, invWbr, rhoGas, ix1, ix2, jy1, jy2, tol)
         implicit none
         integer, intent(in) :: ix1, ix2, jy1, jy2
         real, intent(in) :: dx, dy, invWbr, rhoGas
         real, dimension(:, :, :), intent(in) :: phi, pf
         real, dimension(:, :, :), intent(inout) :: sigx, sigy, crv
         real, intent(in) :: tol
      end subroutine mph_setWeberJumps2d
   end interface

   interface
      subroutine mph_setWeberJumps3d(phi, crv, pf, sigx, sigy, sigz, dx, dy, dz, invWbr, rhoGas, &
                                     ix1, ix2, jy1, jy2, kz1, kz2, tol)
         implicit none
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
         real, intent(in) :: dx, dy, dz, invWbr, rhoGas
         real, dimension(:, :, :), intent(in) :: phi, crv, pf
         real, dimension(:, :, :), intent(inout) :: sigx, sigy, sigz
         real, intent(in) :: tol
      end subroutine mph_setWeberJumps3d
   end interface

   interface
      subroutine mph_init()
         implicit none
      end subroutine mph_init
   end interface

End module mph_interface
