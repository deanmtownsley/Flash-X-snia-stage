!!****if* source/physics/sourceTerms/Heater/Heater_interface
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
!!  Heater_interface()
!!
!! DESCRIPTION
!!  This is an interface specific for heater geometry and specifications
!!
!!***

#include "constants.h"
#include "Simulation.h"

module Heater_interface

   implicit none

   interface
      subroutine Heater_init()
      end subroutine Heater_init
   end interface

   interface
      subroutine Heater_finalize()
      end subroutine Heater_finalize
   end interface

   interface
      subroutine Heater_checkSites(tileDesc, blockCount)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         integer, intent(in) :: blockCount
      end subroutine Heater_checkSites
   end interface

   interface
      subroutine Heater_lsReInit(tileDesc, stime, blockCount)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         real, intent(in) :: stime
         integer, intent(in) :: blockCount
      end subroutine Heater_lsReInit
   end interface

   interface
      subroutine Heater_applyBCToRegion(level, ivar, gridDataStruct, regionData, coordinates, regionSize, &
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
      end subroutine Heater_applyBCToRegion
   end interface

   interface
      subroutine Heater_tagSites(stime)
         real, intent(in) :: stime
      end subroutine Heater_tagSites
   end interface

   interface
      subroutine Heater_mapSitesToProc(initial, gridChanged)
         implicit none
         logical, intent(in), optional :: initial
         logical, intent(in), optional :: gridChanged
      end subroutine Heater_mapSitesToProc
   end interface

   interface
      subroutine Heater_initBlk(xcell, ycell, zcell, ix1, ix2, jy1, jy2, kz1, kz2, temp, phi)
         real, dimension(:, :, :), intent(inout) :: temp
         real, dimension(:, :, :), intent(inout), optional :: phi
         real, dimension(:), intent(in)        :: xcell, ycell, zcell
         integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
      end subroutine Heater_initBlk
   end interface

end module Heater_interface
