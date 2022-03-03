!!  Multiphase_interface
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
!! SYNOPSIS
!!
!!  use Multiphase_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Multiphase module
!! module that defines its public interfaces.
!!
!!***

Module Multiphase_interface

   implicit none

#include "constants.h"
#include "Simulation.h"

   interface
      subroutine Multiphase_init(restart)
         implicit none
         logical, intent(in) :: restart
      end subroutine Multiphase_init
   end interface

   interface
      subroutine Multiphase_advection(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine Multiphase_advection
   end interface

   interface
      subroutine Multiphase_solve(tileDesc, dt)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         real, INTENT(IN) :: dt
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine Multiphase_solve
   end interface

   interface
      subroutine Multiphase_redistance(tileDesc, iteration)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         integer, INTENT(IN) :: iteration
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine Multiphase_redistance
   end interface

   interface
      subroutine Multiphase_indicators()
         implicit none
      end subroutine Multiphase_indicators
   end interface

   interface
      subroutine Multiphase_finalize()
         implicit none
      end subroutine Multiphase_finalize
   end interface

   interface
      subroutine Multiphase_setFluidProps(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine Multiphase_setFluidProps
   end interface

   interface
      subroutine Multiphase_setThermalProps(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine Multiphase_setThermalProps
   end interface

   interface
      subroutine Multiphase_setPressureJumps(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine Multiphase_setPressureJumps
   end interface

   interface
      subroutine Multiphase_getScalarProp(name, value)
         implicit none
         character(len=*), intent(in)  :: name
         real, intent(out)             :: value
      end subroutine Multiphase_getScalarProp
   end interface

   interface
      subroutine Multiphase_thermalForcing(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine Multiphase_thermalForcing
   end interface

   interface
      subroutine Multiphase_divergence(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine Multiphase_divergence
   end interface

   interface
      subroutine Multiphase_extrapFluxes(tileDesc, iteration)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         integer, INTENT(IN) :: iteration
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine Multiphase_extrapFluxes
   end interface

   interface
      subroutine Multiphase_setMassFlux(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine Multiphase_setMassFlux
   end interface

   interface
      subroutine Multiphase_velForcing(tileDesc, dt)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         real, intent(in) :: dt
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine Multiphase_velForcing
   end interface

   interface
      subroutine Multiphase_reInitGridVars(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine Multiphase_reInitGridVars
   end interface

   interface
      subroutine Multiphase_getGridVar(name, value)
         implicit none
         character(len=*), intent(in)  :: name
         integer, intent(out)          :: value
      end subroutine Multiphase_getGridVar
   end interface

end module Multiphase_interface
