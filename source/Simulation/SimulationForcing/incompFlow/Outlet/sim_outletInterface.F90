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
      subroutine sim_outletLSDampingBlk2d(pfrc, phi, xcenter, ycenter, boundBox, &
                                          dt, dx, dy, ix1, ix2, jy1, jy2, &
                                          outletFlag, outletSink, outletBuffer, outletGrowthRate, &
                                          xMin, xMax, yMin, yMax)
         real, dimension(:, :, :), intent(inout) :: pfrc
         real, dimension(:, :, :), intent(in)    :: phi
         real, dimension(:), intent(in)          :: xcenter, ycenter
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy
         integer, intent(in)                     :: ix1, ix2, jy1, jy2
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag
         real, intent(in) :: outletSink, outletBuffer, outletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax
      end subroutine sim_outletLSDampingBlk2d

      subroutine sim_outletLSDampingBlk3d(pfrc, phi, xcenter, ycenter, zcenter, boundBox, &
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
      end subroutine sim_outletLSDampingBlk3d
   end interface

   interface sim_outletVelBlk
      subroutine sim_outletVelBlk2d(u, v, ru, rv, xcenter, ycenter, &
                                    boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, &
                                    inletFlag, inletBuffer, inletGrowthRate, &
                                    outletFlag, outletBuffer, outletGrowthRate, &
                                    volAux, QAux, QOut, &
                                    xMin, xMax, yMin, yMax, gravX, gravY)

         implicit none
         real, dimension(:, :, :), intent(in)    :: u, v
         real, dimension(:, :, :), intent(inout) :: ru, rv
         real, dimension(:), intent(in)          :: xcenter, ycenter
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy
         integer, intent(in)                     :: ix1, ix2, jy1, jy2
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag, inletFlag
         real, intent(in) :: inletBuffer, inletGrowthRate, outletBuffer, outletGrowthRate
         real, intent(inout) :: QAux(LOW:HIGH, MDIM), volAux(LOW:HIGH, MDIM)
         real, intent(in) :: QOut(LOW:HIGH, MDIM)
         real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY
      end subroutine sim_outletVelBlk2d

      subroutine sim_outletVelBlk3d(u, v, w, ru, rv, rw, xcenter, ycenter, zcenter, &
                                    boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                    inletFlag, inletBuffer, inletGrowthRate, &
                                    outletFlag, outletBuffer, outletGrowthRate, &
                                    volAux, QAux, QOut, &
                                    xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ)

         implicit none
         real, dimension(:, :, :), intent(in)  :: u, v, w
         real, dimension(:, :, :), intent(inout) :: ru, rv, rw
         real, dimension(:), intent(in)        :: xcenter, ycenter, zcenter
         real, dimension(:, :), intent(in)     :: boundBox
         real, intent(in)                      :: dt, dx, dy, dz
         integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag, inletFlag
         real, intent(in) :: inletBuffer, inletGrowthRate, outletBuffer, outletGrowthRate
         real, intent(inout) :: QAux(LOW:HIGH, MDIM), volAux(LOW:HIGH, MDIM)
         real, intent(in) :: QOut(LOW:HIGH, MDIM)
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ
      end subroutine sim_outletVelBlk3d
   end interface

   interface sim_outletVelBlkPhased
      subroutine sim_outletVelBlk2dPhased(u, v, ru, rv, phi, xcenter, ycenter, &
                                          boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, &
                                          inletFlag, inletBuffer, inletGrowthRate, &
                                          outletFlag, outletBuffer, outletGrowthRate, &
                                          volAuxLiq, volAuxGas, QAuxLiq, QAuxGas, QOutLiq, QOutGas, &
                                          xMin, xMax, yMin, yMax, gravX, gravY)

         implicit none
         real, dimension(:, :, :), intent(in)    :: u, v
         real, dimension(:, :, :), intent(inout) :: ru, rv
         real, dimension(:, :, :), intent(in) :: phi
         real, dimension(:), intent(in)          :: xcenter, ycenter
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy
         integer, intent(in)                     :: ix1, ix2, jy1, jy2
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag, inletFlag
         real, intent(in) :: inletBuffer, inletGrowthRate, outletBuffer, outletGrowthRate
         real, intent(inout) :: QAuxLiq(LOW:HIGH, MDIM), QAuxGas(LOW:HIGH, MDIM)
         real, intent(inout) :: volAuxLiq(LOW:HIGH, MDIM), volAuxGas(LOW:HIGH, MDIM)
         real, intent(in) :: QOutLiq(LOW:HIGH, MDIM), QOutGas(LOW:HIGH, MDIM)
         real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY
      end subroutine sim_outletVelBlk2dPhased

      subroutine sim_outletVelBlk3dPhased(u, v, w, ru, rv, rw, phi, xcenter, ycenter, zcenter, &
                                          boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                          inletFlag, inletBuffer, inletGrowthRate, &
                                          outletFlag, outletBuffer, outletGrowthRate, &
                                          volAuxLiq, volAuxGas, QAuxLiq, QAuxGas, QOutLiq, QOutGas, &
                                          xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ)

         implicit none
         real, dimension(:, :, :), intent(in)  :: u, v, w
         real, dimension(:, :, :), intent(inout) :: ru, rv, rw
         real, dimension(:, :, :), intent(in) :: phi
         real, dimension(:), intent(in)        :: xcenter, ycenter, zcenter
         real, dimension(:, :), intent(in)     :: boundBox
         real, intent(in)                      :: dt, dx, dy, dz
         integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag, inletFlag
         real, intent(in) :: inletBuffer, inletGrowthRate, outletBuffer, outletGrowthRate
         real, intent(inout) :: QAuxLiq(LOW:HIGH, MDIM), QAuxGas(LOW:HIGH, MDIM)
         real, intent(inout) :: volAuxLiq(LOW:HIGH, MDIM), volAuxGas(LOW:HIGH, MDIM)
         real, intent(in) :: QOutLiq(LOW:HIGH, MDIM), QOutGas(LOW:HIGH, MDIM)
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ
      end subroutine sim_outletVelBlk3dPhased
   end interface

   interface
      subroutine sim_outletApplyBCToFace(level, ivar, gridDataStruct, regionData, coordinates, regionSize, &
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

      end subroutine sim_outletApplyBCToFace
   end interface

End module sim_outletInterface
