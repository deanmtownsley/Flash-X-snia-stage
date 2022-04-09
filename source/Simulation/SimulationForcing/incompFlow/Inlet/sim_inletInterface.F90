!!****if* source/Simulation/SimulationForcing/incompFlow/Inlet/sim_inletInterface
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
!!  sim_inletInterface()
!!
!! DESCRIPTION
!!  This is an interface specific for inlet boundary conditions
!!
!!***

#include "constants.h"
#include "Simulation.h"

Module sim_inletInterface

   implicit none

   interface
      subroutine sim_inletInit()
      end subroutine sim_inletInit
   end interface

   interface
      subroutine sim_inletFinalize()
      end subroutine sim_inletFinalize
   end interface

   interface
      subroutine sim_inletSetForcing(tileDesc, dt)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         real, intent(in) :: dt
      end subroutine sim_inletSetForcing
   end interface

   interface
      subroutine sim_inletVelBlk2d(u, v, ru, rv, xcell, ycell, &
                                   boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, &
                                   inletFlag, outletFlag, inletBuffer, inletGrowthRate, &
                                   xMin, xMax, yMin, yMax, gravX, gravY)

         implicit none
         real, dimension(:, :, :), intent(in)    :: u, v
         real, dimension(:, :, :), intent(inout) :: ru, rv
         real, dimension(:), intent(in)          :: xcell, ycell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy
         integer, intent(in)                     :: ix1, ix2, jy1, jy2
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: inletFlag, outletFlag
         real, intent(in) :: inletBuffer, inletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY

      end subroutine sim_inletVelBlk2d

      subroutine sim_inletVelBlk3d(u, v, w, ru, rv, rw, xcell, ycell, zcell, &
                                   boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                   inletFlag, outletFlag, inletBuffer, inletGrowthRate, &
                                   xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ)

         implicit none
         real, dimension(:, :, :), intent(in)  :: u, v, w
         real, dimension(:, :, :), intent(inout) :: ru, rv, rw
         real, dimension(:), intent(in)        :: xcell, ycell, zcell
         real, dimension(:, :), intent(in)     :: boundBox
         real, intent(in)                      :: dt, dx, dy, dz
         integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: inletFlag, outletFlag
         real, intent(in) :: inletBuffer, inletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ

      end subroutine sim_inletVelBlk3d
   end interface

   interface
      subroutine sim_inletVelBlk2dPhased(u, v, ru, rv, phi, xcell, ycell, &
                                         boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, &
                                         inletFlag, outletFlag, inletBuffer, inletGrowthRate, &
                                         xMin, xMax, yMin, yMax, gravX, gravY)

         implicit none
         real, dimension(:, :, :), intent(in)    :: u, v
         real, dimension(:, :, :), intent(inout) :: ru, rv
         real, dimension(:, :, :), intent(in) :: phi
         real, dimension(:), intent(in)          :: xcell, ycell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy
         integer, intent(in)                     :: ix1, ix2, jy1, jy2
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag, inletFlag
         real, intent(in) :: inletBuffer, inletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY

      end subroutine sim_inletVelBlk2dPhased

      subroutine sim_inletVelBlk3dPhased(u, v, w, ru, rv, rw, phi, xcell, ycell, zcell, &
                                         boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                         inletFlag, outletFlag, inletBuffer, inletGrowthRate, &
                                         xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ)

         implicit none
         real, dimension(:, :, :), intent(in)  :: u, v, w
         real, dimension(:, :, :), intent(inout) :: ru, rv, rw
         real, dimension(:, :, :), intent(in) :: phi
         real, dimension(:), intent(in)        :: xcell, ycell, zcell
         real, dimension(:, :), intent(in)     :: boundBox
         real, intent(in)                      :: dt, dx, dy, dz
         integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: outletFlag, inletFlag
         real, intent(in) :: inletBuffer, inletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ

      end subroutine sim_inletVelBlk3dPhased
   end interface

   interface
      subroutine sim_inletApplyBCToFace(level, ivar, gridDataStruct, regionData, coordinates, regionSize, &
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

      end subroutine sim_inletApplyBCToFace
   end interface

   interface
      subroutine sim_inletLSDampingBlk2d(pfrc, phi, xcell, ycell, boundBox, &
                                         dt, dx, dy, ix1, ix2, jy1, jy2, &
                                         inletFlag, inletSink, inletBuffer, inletGrowthRate, &
                                         xMin, xMax, yMin, yMax)

         implicit none

         real, dimension(:, :, :), intent(inout) :: pfrc
         real, dimension(:, :, :), intent(in) :: phi
         real, dimension(:), intent(in) :: xcell, ycell
         real, dimension(:, :), intent(in) :: boundBox
         real, intent(in) :: dt, dx, dy
         integer, intent(in) :: ix1, ix2, jy1, jy2
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: inletFlag
         real, intent(in) :: inletSink, inletBuffer, inletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax

      end subroutine sim_inletLSDampingBlk2d

      subroutine sim_inletLSDampingBlk3d(pfrc, phi, xcell, ycell, zcell, boundBox, &
                                         dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                         inletFlag, inletSink, inletBuffer, inletGrowthRate, &
                                         xMin, xMax, yMin, yMax, zMin, zMax)

         implicit none

         real, dimension(:, :, :), intent(inout) :: pfrc
         real, dimension(:, :, :), intent(in) :: phi
         real, dimension(:), intent(in) :: xcell, ycell, zcell
         real, dimension(:, :), intent(in) :: boundBox
         real, intent(in) :: dt, dx, dy, dz
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
         integer, dimension(LOW:HIGH, MDIM), intent(in) :: inletFlag
         real, intent(in) :: inletSink, inletBuffer, inletGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax

      end subroutine sim_inletLSDampingBlk3d
   end interface

End module sim_inletInterface
