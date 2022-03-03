!!****if* source/Simulation/SimulationMain/incompFlow/PoolBoiling/localAPI/heater/sim_heaterInterface
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
      subroutine sim_heaterCheckSites()
      end subroutine sim_heaterCheckSites
   end interface

   interface
      subroutine sim_heaterLSReInit(stime)
         real, intent(in) :: stime
      end subroutine sim_heaterLSReInit
   end interface

   interface
      subroutine sim_heaterApplyBC(dt)
         real, intent(in) :: dt
      end subroutine sim_heaterApplyBC
   end interface

   interface
      subroutine sim_heaterRead(heaterID, heaterFile)
         integer, intent(in)          :: heaterID
         character(len=*), intent(in) :: heaterFile
      end subroutine sim_heaterRead
   end interface

   interface
      subroutine sim_heaterInitBlk(phi, temp, xcell, ycell, zcell, ix1, ix2, jy1, jy2, kz1, kz2)
         real, dimension(:, :, :), intent(inout) :: phi, temp
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

   interface sim_heaterApplyBCToBlk
      subroutine sim_heaterApplyBCToBlk2d(pfrc, tfrc, phi, temp, xcell, ycell, boundBox, dt, dx, dy, ix1, ix2, jy1, jy2)
         real, dimension(:, :, :), intent(inout)  :: pfrc, tfrc
         real, dimension(:, :, :), intent(in)     :: phi, temp
         real, dimension(:), intent(in)         :: xcell, ycell
         real, dimension(:, :), intent(in)       :: boundBox
         real, intent(in)                       :: dt, dx, dy
         integer, intent(in)                    :: ix1, ix2, jy1, jy2
      end subroutine sim_heaterApplyBCToBlk2d
   end interface

   interface
      subroutine sim_heaterApplyBCToFace(level, ivar, gridDataStruct, regionData, coordinates, regionSize, &
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
