!!****if* source/Simulation/SimulationForcing/incompFlow/Heater/Heater_Interface
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
!!  Heater_Interface()
!!
!! DESCRIPTION
!!  This is an interface specific for heater geometry and specifications
!!
!!***

#include "constants.h"
#include "Simulation.h"

Module Heater_interface

   implicit none

   interface
      subroutine Heater_Init()
      end subroutine Heater_Init
   end interface

   interface
      subroutine Heater_Finalize()
      end subroutine Heater_Finalize
   end interface

   interface
      subroutine Heater_CheckSites(tileDesc, blockCount)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         integer, intent(in) :: blockCount
      end subroutine Heater_CheckSites
   end interface

   interface
      subroutine Heater_LSReInit(tileDesc, stime, blockCount)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         real, intent(in) :: stime
         integer, intent(in) :: blockCount
      end subroutine Heater_LSReInit
   end interface

   interface
      subroutine Heater_Read(heaterID, heaterFile)
         integer, intent(in)          :: heaterID
         character(len=*), intent(in) :: heaterFile
      end subroutine Heater_Read
   end interface

   interface
      subroutine Heater_InitBlk(xcell, ycell, zcell, ix1, ix2, jy1, jy2, kz1, kz2, temp, phi)
         real, dimension(:, :, :), intent(inout) :: temp
         real, dimension(:, :, :), intent(inout), optional :: phi
         real, dimension(:), intent(in)        :: xcell, ycell, zcell
         integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
      end subroutine Heater_InitBlk
   end interface

   interface
      subroutine Heater_LSReInitBlk(phi, xcell, ycell, zcell, boundBox, stime, &
                                       ix1, ix2, jy1, jy2, kz1, kz2, lblock)
         real, dimension(:, :, :), intent(inout) :: phi
         real, dimension(:), intent(in)        :: xcell, ycell, zcell
         real, dimension(:, :), intent(in)      :: boundBox
         real, intent(in)                      :: stime
         integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2, lblock
      end subroutine Heater_LSReInitBlk
   end interface

   interface Heater_CheckSitesBlk
      subroutine Heater_CheckSitesBlk2d(phi, xcell, ycell, boundBox, ix1, ix2, jy1, jy2, lblock)
         real, dimension(:, :, :), intent(in)  :: phi
         real, dimension(:), intent(in)      :: xcell, ycell
         real, dimension(:, :), intent(in)    :: boundBox
         integer, intent(in)                 :: ix1, ix2, jy1, jy2, lblock
      end subroutine Heater_CheckSitesBlk2d

      subroutine Heater_CheckSitesBlk3d(phi, xcell, ycell, zcell, boundBox, ix1, ix2, jy1, jy2, kz1, kz2, lblock)
         real, dimension(:, :, :), intent(in)  :: phi
         real, dimension(:), intent(in)      :: xcell, ycell, zcell
         real, dimension(:, :), intent(in)    :: boundBox
         integer, intent(in)                :: ix1, ix2, jy1, jy2, kz1, kz2, lblock
      end subroutine Heater_CheckSitesBlk3d
   end interface

   interface
      subroutine Heater_ApplyBCToRegion(level, ivar, gridDataStruct, regionData, coordinates, regionSize, &
                                           guard, face, axis, secondDir, thirdDir)
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
      end subroutine
   end interface

   interface
      subroutine Heater_TagSites(stime)
         real, intent(in) :: stime
      end subroutine
   end interface

   interface
      subroutine sim_forceHeater(nstep, dt, stime)
         implicit none
         integer, intent(in) :: nstep
         real, intent(in) :: dt
         real, intent(in) :: stime
      end subroutine sim_forceHeater
   end interface

   interface
      subroutine Heater_MapSitesToProc(initial, gridChanged)
         implicit none
         logical, intent(in), optional :: initial
         logical, intent(in), optional :: gridChanged
      end subroutine Heater_MapSitesToProc
   end interface

#ifdef SIM_HEATER_ANN_SEARCH
   interface
      subroutine Heater_AnnBuildTree(heater)
         use Heater_Data, ONLY: Heater_Type
         type(Heater_Type), intent(INOUT)  :: heater
      end subroutine Heater_AnnBuildTree
   end interface

   interface
      subroutine Heater_AnnSearchTree(heater, queryPt, annElems, annIdx)
         use Heater_Data, ONLY: Heater_Type
         type(Heater_Type), intent(IN)  :: heater
         integer, intent(IN) :: annElems
         ! query point
         real, dimension(:), target, intent(IN) :: queryPt
         ! indices of nearest neighbors
         integer, dimension(:), target, intent(OUT):: annIdx
      end subroutine Heater_AnnSearchTree
   end interface
#endif

End module Heater_interface
