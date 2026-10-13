!****h* source/physics/ImBound/Imbound_interface
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
!!  Imbound_interface
!!
!! SYNOPSIS
!!
!!  use Imbound_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Immersed boundary (IB)
!! unit that defines its public interfaces.
!!
!!***

Module ImBound_interface

   implicit none

#include "Simulation.h"

   interface !ImBound_init
      subroutine ImBound_init(restart)
         implicit none
         logical, INTENT(IN) :: restart
      end subroutine ImBound_init
   end interface

   interface  !ImBound_finalize
      subroutine ImBound_finalize()
         implicit none
      end subroutine ImBound_finalize
   end interface

   interface
      subroutine ImBound_getScalarProp(name, value)
         implicit none
         character(len=*), intent(in)  :: name
         real, intent(out)             :: value
      end subroutine ImBound_getScalarProp
   end interface

   interface
      subroutine ImBound_mapToGrid(tileDesc, bodyInfo)
         use Grid_tile, ONLY: Grid_tile_t
         use ImBound_type, ONLY: ImBound_type_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
         type(ImBound_type_t), intent(in) :: bodyInfo
      end subroutine ImBound_mapToGrid
   end interface

   interface
      subroutine ImBound_getBodyPtr(bodyPtr, ibd)
         use ImBound_type, only: ImBound_type_t
         type(ImBound_type_t), pointer :: bodyPtr
         integer, intent(in) :: ibd
      end subroutine ImBound_getBodyPtr
   end interface

   interface
      subroutine ImBound_releaseBodyPtr(bodyPtr, ibd)
         use ImBound_type, only: ImBound_type_t
         type(ImBound_type_t), pointer :: bodyPtr
         integer, intent(in) :: ibd
      end subroutine ImBound_releaseBodyPtr
   end interface

   interface
      subroutine ImBound_redistanceLS(tileDesc, iteration)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         integer, intent(in) :: iteration
      end subroutine ImBound_redistanceLS
   end interface

   interface
      subroutine ImBound_advectLS(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
      end subroutine ImBound_advectLS
   end interface

   interface
      subroutine ImBound_solveLS(tileDesc, dt)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         real, intent(in) :: dt
      end subroutine ImBound_solveLS
   end interface

   interface
      subroutine ImBound_reInitGridVars(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine ImBound_reInitGridVars
   end interface

   interface
      subroutine ImBound_advance(bodyInfo, time, dt)
         use ImBound_type, ONLY: ImBound_type_t
         implicit none
         type(ImBound_type_t), intent(inout) :: bodyInfo
         real, intent(in) :: time, dt
      end subroutine ImBound_advance
   end interface

   interface
      subroutine ImBound_skipBox(tileDesc, bodyInfo, skipBox)
         use Grid_tile, ONLY: Grid_tile_t
         use ImBound_type, ONLY: ImBound_type_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
         type(ImBound_type_t), intent(in) :: bodyInfo
         logical, intent(out) :: skipBox
      end subroutine ImBound_skipBox
   end interface
   
   interface ImBound_velForcing
      subroutine ImBound_velForcing_fixed(tileDesc, dt)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         real, intent(in) :: dt
         type(Grid_tile_t), intent(in) :: tileDesc
      end subroutine ImBound_velForcing_fixed

      subroutine ImBound_velForcing_moving(tileDesc, bodyInfo, dt)
         use Grid_tile, ONLY: Grid_tile_t
         use ImBound_type, ONLY: ImBound_type_t
         implicit none
         real, intent(in) :: dt
         type(Grid_tile_t), intent(in) :: tileDesc
         type(ImBound_type_t), intent(in) :: bodyInfo
      end subroutine ImBound_velForcing_moving
   end interface ImBound_velForcing

   interface
      subroutine ImBound_initBlk(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
      end subroutine ImBound_initBlk
   end interface

end Module ImBound_interface
