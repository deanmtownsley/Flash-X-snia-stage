!!****if* source/Simulation/SimulationForcing/incompFlow/Outlet/sim_outletInterface
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
!! NAME
!!
!!
!! SYNOPSIS
!!  sim_outletInterface()
!!
!! DESCRIPTION
!!  This is an interface specific for outlet boundary conditions
!!
!!***

#include "constants.h"
#include "Simulation.h"

Module sim_outletInterface

   implicit none

   interface
      subroutine sim_outletInit()
      end subroutine sim_outletInit
   end interface

   interface
      subroutine sim_outletFinalize()
      end subroutine sim_outletFinalize
   end interface

   interface
      subroutine sim_outletSetForcing(tileDesc, dt)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         real, intent(in) :: dt
      end subroutine sim_outletSetForcing
   end interface

   interface
      subroutine sim_outletLSDamping(pfrc, phi, xcenter, ycenter, zcenter, boundBox, &
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
      end subroutine sim_outletLSDamping
   end interface

   interface
      subroutine sim_outletVelFrc(vel, rhs, xgrid, ygrid, zgrid, &
                                  dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                  xMin, xMax, yMin, yMax, zMin, zMax, &
                                  outletFlag, outletBuffer, outletGrowthRate, &
                                  axis, volAux, QAux, QOut)
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
      end subroutine sim_outletVelFrc
   end interface

   interface
      subroutine sim_outletVelFrcPhased(vel, rhs, phi, xgrid, ygrid, zgrid, &
                                        dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                        xMin, xMax, yMin, yMax, zMin, zMax, &
                                        outletFlag, outletBuffer, outletGrowthRate, &
                                        axis, volAuxLiq, volAuxGas, QAuxLiq, QAuxGas, &
                                        QOutLiq, QOutGas)

         implicit none
         real, dimension(:, :, :), intent(in) :: vel, phi
         real, dimension(:, :, :), intent(inout) :: rhs
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
      end subroutine sim_outletVelFrcPhased
   end interface

   interface
      subroutine sim_outletApplyBCToRegion(level, ivar, gridDataStruct, regionData, coordinates, regionSize, &
                                         guard, face, axis, secondDir, thirdDir)

         implicit none
         integer, intent(IN) :: level, ivar, gridDataStruct
         integer, dimension(REGION_DIM), intent(IN) :: regionSize
         real, dimension(regionSize(BC_DIR), &
                         regionSize(SECOND_DIR), &
                         regionSize(THIRD_DIR), &
                         regionSize(STRUCTSIZE)), intent(INOUT) :: regionData
         real, dimension(regionSize(BC_DIR), &
                         regionSize(SECOND_DIR), &
                         regionSize(THIRD_DIR), &
                         MDIM), intent(IN) :: coordinates
         integer, intent(IN) :: guard, face, axis, secondDir, thirdDir

      end subroutine sim_outletApplyBCToRegion
   end interface

End module sim_outletInterface
