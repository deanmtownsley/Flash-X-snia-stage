!!****if* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterInterface
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
!!  sim_heaterInterface()
!!
!! DESCRIPTION
!!  This is an interface specific for heater geometry and specifications
!!
!!***

#include "constants.h"

Module sim_heaterInterface

   implicit none

   interface
      subroutine sim_heaterInit()
      end subroutine sim_heaterInit
   end interface

   interface
      subroutine sim_heaterFinalize()
      end subroutine sim_heaterFinalize
   end interface

   interface
      subroutine sim_heaterCheckSites(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
      end subroutine sim_heaterCheckSites
   end interface

   interface
      subroutine sim_heaterLSReInit(tileDesc, stime)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         real, intent(in) :: stime
      end subroutine sim_heaterLSReInit
   end interface

   interface
      subroutine sim_heaterRead(heaterID, heaterFile)
         integer, intent(in)          :: heaterID
         character(len=*), intent(in) :: heaterFile
      end subroutine sim_heaterRead
   end interface

   interface
      subroutine sim_heaterInitBlk(xcell, ycell, zcell, ix1, ix2, jy1, jy2, kz1, kz2, temp, phi)
         real, dimension(:, :, :), intent(inout) :: temp
         real, dimension(:, :, :), intent(inout), optional :: phi
         real, dimension(:), intent(in)        :: xcell, ycell, zcell
         integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
      end subroutine sim_heaterInitBlk
   end interface

   interface
      subroutine sim_heaterLSReInitBlk(phi, xcell, ycell, zcell, boundBox, stime, ix1, ix2, jy1, jy2, kz1, kz2)
         real, dimension(:, :, :), intent(inout) :: phi
         real, dimension(:), intent(in)        :: xcell, ycell, zcell
         real, dimension(:, :), intent(in)      :: boundBox
         real, intent(in)                      :: stime
         integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
      end subroutine sim_heaterLSReInitBlk
   end interface

   interface sim_heaterCheckSitesBlk
      subroutine sim_heaterCheckSitesBlk2d(phi, xcell, ycell, boundBox, ix1, ix2, jy1, jy2)
         real, dimension(:, :, :), intent(in)  :: phi
         real, dimension(:), intent(in)      :: xcell, ycell
         real, dimension(:, :), intent(in)    :: boundBox
         integer, intent(in)                 :: ix1, ix2, jy1, jy2
      end subroutine sim_heaterCheckSitesBlk2d

      subroutine sim_heaterCheckSitesBlk3d(phi, xcell, ycell, zcell, boundBox, ix1, ix2, jy1, jy2, kz1, kz2)
         real, dimension(:, :, :), intent(in)  :: phi
         real, dimension(:), intent(in)      :: xcell, ycell, zcell
         real, dimension(:, :), intent(in)    :: boundBox
         integer, intent(in)                :: ix1, ix2, jy1, jy2, kz1, kz2
      end subroutine sim_heaterCheckSitesBlk3d
   end interface

   interface
      subroutine sim_heaterApplyBCToRegion(level, ivar, gridDataStruct, regionData, coordinates, regionSize, &
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
      subroutine sim_heaterTagSites(stime)
         real, intent(in) :: stime
      end subroutine
   end interface

End module sim_heaterInterface
