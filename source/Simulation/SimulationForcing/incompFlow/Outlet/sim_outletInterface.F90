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

   interface sim_outletLSDampingBlk
      subroutine sim_outletLSDampingBlk2d(pfrc, phi, xcell, ycell, boundBox, &
                                          dt, dx, dy, ix1, ix2, jy1, jy2, &
                                          outletFlag, outletSink, outletBuffer, outletGrowthRate, &
                                          xMin, xMax, yMin, yMax)
         real, dimension(:, :, :), intent(inout) :: pfrc
         real, dimension(:, :, :), intent(in)    :: phi
         real, dimension(:), intent(in)          :: xcell, ycell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy
         integer, intent(in)                     :: ix1, ix2, jy1, jy2
         integer, dimension(2, MDIM), intent(in) :: outletFlag
         real, intent(in) :: outletSink, outletBuffer, outletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax
      end subroutine sim_outletLSDampingBlk2d

      subroutine sim_outletLSDampingBlk3d(pfrc, phi, xcell, ycell, zcell, boundBox, &
                                          dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                          outletFlag, outletSink, outletBuffer, outletGrowthRate, &
                                          xMin, xMax, yMin, yMax, zMin, zMax)
         real, dimension(:, :, :), intent(inout) :: pfrc
         real, dimension(:, :, :), intent(in)    :: phi
         real, dimension(:), intent(in)          :: xcell, ycell, zcell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy, dz
         integer, intent(in)                     :: ix1, ix2, jy1, jy2, kz1, kz2
         integer, dimension(2, MDIM), intent(in) :: outletFlag
         real, intent(in) :: outletSink, outletBuffer, outletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax
      end subroutine sim_outletLSDampingBlk3d
   end interface

   interface sim_outletVelBlk
      subroutine sim_outletVelBlk2d(u, v, ru, rv, xcell, ycell, &
                                    boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, inletFlag, &
                                    outletFlag, velAux, velOut, outletBuffer, outletGrowthRate, &
                                    xMin, xMax, yMin, yMax, gravX, gravY)

         real, dimension(:, :, :), intent(in)    :: u, v
         real, dimension(:, :, :), intent(inout) :: ru, rv
         real, dimension(:), intent(in)          :: xcell, ycell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy
         integer, intent(in)                     :: ix1, ix2, jy1, jy2
         integer, dimension(2, MDIM), intent(in) :: outletFlag, inletFlag
         real, intent(inout) :: velAux(2, MDIM)
         real, intent(in) :: velOut(2, MDIM), outletBuffer, outletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY
      end subroutine sim_outletVelBlk2d

      subroutine sim_outletVelBlk3d(u, v, w, ru, rv, rw, xcell, ycell, zcell, &
                                    boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, inletFlag, &
                                    outletFlag, velAux, velOut, outletBuffer, outletGrowthRate, &
                                    xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ)

         real, dimension(:, :, :), intent(in)    :: u, v, w
         real, dimension(:, :, :), intent(inout) :: ru, rv, rw
         real, dimension(:), intent(in)          :: xcell, ycell, zcell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy, dz
         integer, intent(in)                     :: ix1, ix2, jy1, jy2, kz1, kz2
         integer, dimension(2, MDIM), intent(in) :: outletFlag, inletFlag
         real, intent(inout) :: velAux(2, MDIM)
         real, intent(in) :: velOut(2, MDIM), outletBuffer, outletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ
      end subroutine sim_outletVelBlk3d
   end interface

   interface sim_outletVelBlkPhased
      subroutine sim_outletVelBlk2dPhased(u, v, ru, rv, phi, xcell, ycell, &
                                          boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, inletFlag, &
                                          outletFlag, phaseAuxLiq, phaseAuxGas, &
                                          velAuxLiq, velAuxGas, velOutLiq, velOutGas, &
                                          outletBuffer, outletGrowthRate, &
                                          xMin, xMax, yMin, yMax, gravX, gravY)

         real, dimension(:, :, :), intent(in)    :: u, v
         real, dimension(:, :, :), intent(inout) :: ru, rv
         real, dimension(:, :, :), intent(in) :: phi
         real, dimension(:), intent(in)          :: xcell, ycell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy
         integer, intent(in)                     :: ix1, ix2, jy1, jy2
         integer, dimension(2, MDIM), intent(in) :: outletFlag, inletFlag
         real, intent(inout) :: velAuxLiq(2, MDIM), velAuxGas(2, MDIM), phaseAuxLiq(2, MDIM), phaseAuxGas(2, MDIM)
         real, intent(in) :: velOutLiq(2, MDIM), velOutGas(2, MDIM), outletBuffer, outletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY
      end subroutine sim_outletVelBlk2dPhased

      subroutine sim_outletVelBlk3dPhased(u, v, w, ru, rv, rw, phi, xcell, ycell, zcell, &
                                          boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, inletFlag, &
                                          outletFlag, phaseAuxLiq, phaseAuxGas, &
                                          velAuxLiq, velAuxGas, velOutLiq, velOutGas, &
                                          outletBuffer, outletGrowthRate, &
                                          xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ)

         real, dimension(:, :, :), intent(in)    :: u, v, w
         real, dimension(:, :, :), intent(inout) :: ru, rv, rw
         real, dimension(:, :, :), intent(in) :: phi
         real, dimension(:), intent(in)          :: xcell, ycell, zcell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy, dz
         integer, intent(in)                     :: ix1, ix2, jy1, jy2, kz1, kz2
         integer, dimension(2, MDIM), intent(in) :: outletFlag, inletFlag
         real, intent(inout) :: velAuxLiq(2, MDIM), velAuxGas(2, MDIM), phaseAuxLiq(2, MDIM), phaseAuxGas(2, MDIM)
         real, intent(in) :: velOutLiq(2, MDIM), velOutGas(2, MDIM), outletBuffer, outletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ
      end subroutine sim_outletVelBlk3dPhased
   end interface

End module sim_outletInterface
