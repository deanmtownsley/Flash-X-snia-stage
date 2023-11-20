!!****if* source/physics/sourceTerms/Outlet/localAPI/out_interface
!!
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  out_interface()
!!
!! DESCRIPTION
!!  This is an interface specific for outlet boundary conditions
!!
!!***

#include "constants.h"
#include "Simulation.h"

Module out_interface

   implicit none

   interface
      subroutine out_lsDamping(pfrc, phi, xcenter, ycenter, zcenter, boundBox, &
                               dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                               outletFlag, outletSink, outletBuffer, outletGrowthRate, &
                               xMin, xMax, yMin, yMax, zMin, zMax)
         real, dimension(:, :, :), intent(inout) :: pfrc
         real, dimension(:, :, :), intent(in)    :: phi
         real, dimension(:), intent(in)          :: xcenter, ycenter, zcenter
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy, dz
         integer, intent(in)                     :: ix1, ix2, jy1, jy2, kz1, kz2
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag
         real, intent(in) :: outletSink, outletBuffer, outletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax
      end subroutine out_lsDamping
   end interface

   interface
      subroutine out_velFrc(vel, rhs, xgrid, ygrid, zgrid, &
                            dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                            xMin, xMax, yMin, yMax, zMin, zMax, &
                            outletFlag, outletBuffer, outletGrowthRate, &
                            axis, volAux, QAux, QOut, velref)
         implicit none
         real, dimension(:, :, :), intent(in) :: vel
         real, dimension(:, :, :), intent(inout) :: rhs
         real, dimension(:), intent(in) :: xgrid, ygrid, zgrid
         real, intent(in) :: dt, dx, dy, dz
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag
         real, intent(in) :: outletBuffer, outletGrowthRate
         integer, intent(in) :: axis
         real, intent(inout) :: QAux(LOW:HIGH, MDIM), volAux(LOW:HIGH, MDIM)
         real, intent(in) :: QOut(LOW:HIGH, MDIM)
         real, intent(in) :: velref
      end subroutine out_velFrc
   end interface

   interface
      subroutine out_velFrcPhased(vel, rhs, sigm, phi, xgrid, ygrid, zgrid, &
                                  dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                  xMin, xMax, yMin, yMax, zMin, zMax, &
                                  outletFlag, outletBuffer, outletGrowthRate, &
                                  axis, volAuxLiq, volAuxGas, QAuxLiq, QAuxGas, &
                                  QOutLiq, QOutGas, velref)

         implicit none
         real, dimension(:, :, :), intent(in) :: vel, phi
         real, dimension(:, :, :), intent(inout) :: rhs, sigm
         real, dimension(:), intent(in) :: xgrid, ygrid, zgrid
         real, intent(in) :: dt, dx, dy, dz
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag
         real, intent(in) :: outletBuffer, outletGrowthRate
         integer, intent(in) :: axis
         real, intent(inout) :: QAuxLiq(LOW:HIGH, MDIM), QAuxGas(LOW:HIGH, MDIM)
         real, intent(inout) :: volAuxLiq(LOW:HIGH, MDIM), volAuxGas(LOW:HIGH, MDIM)
         real, intent(in) :: QOutLiq(LOW:HIGH, MDIM), QOutGas(LOW:HIGH, MDIM)
         real, intent(in) :: velref
      end subroutine out_velFrcPhased
   end interface

End module out_interface
