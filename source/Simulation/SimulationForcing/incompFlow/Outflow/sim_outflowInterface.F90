!!****if* source/Simulation/SimulationForcing/incompFlow/Outflow/sim_outflowInterface
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
!!  sim_outflowInterface()
!!
!! DESCRIPTION
!!  This is an interface specific for outflow boundary conditions
!!
!!***

#include "constants.h"
#include "Simulation.h"

Module sim_outflowInterface

   implicit none

   interface
      subroutine sim_outflowInit()
      end subroutine sim_outflowInit
   end interface

   interface
      subroutine sim_outflowFinalize()
      end subroutine sim_outflowFinalize
   end interface

   interface
      subroutine sim_outflowSetForcing(tileDesc, velOutAux, dt)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         real, intent(inout) :: velOutAux(LOW:HIGH, MDIM)
         real, intent(in) :: dt
      end subroutine sim_outflowSetForcing
   end interface

   interface sim_outflowLSDampingBlk
      subroutine sim_outflowLSDampingBlk2d(pfrc, phi, xcell, ycell, boundBox, &
                                           dt, dx, dy, ix1, ix2, jy1, jy2, &
                                           domainBC, outflowSink, outflowBuffer, outflowGrowthRate, &
                                           xMin, xMax, yMin, yMax)
         real, dimension(:, :, :), intent(inout) :: pfrc
         real, dimension(:, :, :), intent(in)    :: phi
         real, dimension(:), intent(in)          :: xcell, ycell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy
         integer, intent(in)                     :: ix1, ix2, jy1, jy2
         integer, dimension(2, MDIM), intent(in) :: domainBC
         real, intent(in) :: outflowSink, outflowBuffer, outflowGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax
      end subroutine sim_outflowLSDampingBlk2d

      subroutine sim_outflowLSDampingBlk3d(pfrc, phi, xcell, ycell, zcell, boundBox, &
                                           dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                           domainBC, outflowSink, outflowBuffer, outflowGrowthRate, &
                                           xMin, xMax, yMin, yMax, zMin, zMax)
         real, dimension(:, :, :), intent(inout) :: pfrc
         real, dimension(:, :, :), intent(in)    :: phi
         real, dimension(:), intent(in)          :: xcell, ycell, zcell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy, dz
         integer, intent(in)                     :: ix1, ix2, jy1, jy2, kz1, kz2
         integer, dimension(2, MDIM), intent(in) :: domainBC
         real, intent(in) :: outflowSink, outflowBuffer, outflowGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax
      end subroutine sim_outflowLSDampingBlk3d
   end interface

   interface sim_outflowVelBlk
      subroutine sim_outflowVelBlk2d(u, v, ru, rv, xcell, ycell, &
                                     boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, &
                                     domainBC, velOutAux, velOut, outflowBuffer, outflowGrowthRate,&
                                     xMin, xMax, yMin, yMax)

         real, dimension(:, :, :), intent(in)    :: u, v
         real, dimension(:, :, :), intent(inout) :: ru, rv
         real, dimension(:), intent(in)          :: xcell, ycell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy
         integer, intent(in)                     :: ix1, ix2, jy1, jy2
         integer, dimension(2, MDIM), intent(in) :: domainBC
         real, intent(inout) :: velOutAux(2,MDIM)
         real, intent(in) :: velOut(2,MDIM), outflowBuffer, outflowGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax
      end subroutine sim_outflowVelBlk2d

      subroutine sim_outflowVelBlk3d(u, v, w, ru, rv, rw, xcell, ycell, zcell, &
                                     boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                                     domainBC, velOutAux, velOut, outflowBuffer, outflowGrowthRate,&
                                     xMin, xMax, yMin, yMax, zMin, zMax)

         real, dimension(:, :, :), intent(in)    :: u, v, w
         real, dimension(:, :, :), intent(inout) :: ru, rv, rw
         real, dimension(:), intent(in)          :: xcell, ycell, zcell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                        :: dt, dx, dy, dz
         integer, intent(in)                     :: ix1, ix2, jy1, jy2, kz1, kz2
         integer, dimension(2, MDIM), intent(in) :: domainBC
         real, intent(inout) :: velOutAux(2,MDIM)
         real, intent(in) :: velOut(2,MDIM), outflowBuffer, outflowGrowthRate
         real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax
      end subroutine sim_outflowVelBlk3d
   end interface

End module sim_outflowInterface
