!!****if* source/physics/sourceTerms/Outlet/Outlet_interface
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
!!  Outlet_interface()
!!
!! DESCRIPTION
!!  This is an interface specific for outlet boundary conditions
!!
!!***

#include "constants.h"
#include "Simulation.h"

Module Outlet_interface

   implicit none

   interface
      subroutine Outlet_init()
      end subroutine Outlet_init
   end interface

   interface
      subroutine Outlet_finalize()
      end subroutine Outlet_finalize
   end interface

   interface
      subroutine Outlet_setForcing(tileDesc, dt)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         real, intent(in) :: dt
      end subroutine Outlet_setForcing
   end interface

   interface
      subroutine Outlet_applyBCToRegion(level, ivar, gridDataStruct, regionData, coordinates, regionSize, &
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

      end subroutine Outlet_applyBCToRegion
   end interface

   interface
      subroutine Outlet_consolidate()
         implicit none
      end subroutine Outlet_consolidate
   end interface

End module Outlet_interface
