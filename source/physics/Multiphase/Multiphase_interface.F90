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
      !! solnData :: source=grid_data, structure_index=[center, 1], RW = [DFUN_VAR], R = [NRMX_VAR, NRMY_VAR]
      !! facexData :: source=grid_data, structure_index=[facex, 1], RW = [VELC_FACE_VAR]
      !! del :: source=tile_deltas
      !! blkLimits :: source=tile_interior
      !! blkLimitsGC :: source=tile_arrayBounds
      !!subroutine Multiphase_advection(solnData, facexData, del, blkLimits, blkLimitsGC)
      subroutine Multiphase_advection(solnData, facexData, faceyData, facezData,del, lo, hi)
        implicit none
        real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
        real,dimension(MDIM), intent(IN) :: del
        integer, dimension(MDIM),intent(IN) :: lo, hi
      end subroutine Multiphase_advection
   end interface

   interface
      subroutine Multiphase_solve(solnData, lo, hi, dt)
        
        real, dimension(:,:,:,:), pointer :: solnData
        real, INTENT(IN) :: dt
        integer, dimension(MDIM), intent(IN) :: lo, hi
        !-----------------------------------------------------------------------------------------
      end subroutine Multiphase_solve
   end interface

   interface
      subroutine Multiphase_redistance(solnData,del,lo,hi, iteration)
        real, pointer, dimension(:, :, :, :) :: solnData
        integer, intent(in) :: iteration
        integer,dimension(MDIM), intent(IN) :: lo, hi
        real,dimension(MDIM),intent(IN) :: del
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
      subroutine Multiphase_setFluidProps(solnData, facexData, faceyData, facezData, del,&
            logc, higc)
        real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
        real,dimension(MDIM), intent(in) :: del
        integer, intent(IN), dimension(MDIM) :: logc, higc
      end subroutine Multiphase_setFluidProps
   end interface

   interface
      subroutine Multiphase_setThermalProps(solnData)
        real,dimension(:,:,:,:), pointer :: solnData
      end subroutine Multiphase_setThermalProps
   end interface

   interface
      subroutine Multiphase_setPressureJumps(solnData, facexData, faceyData, facezData, del, lo, hi)
        implicit none
        integer, dimension(MDIM),intent(IN) :: lo, hi
        real,dimension(MDIM) ::  del
        real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
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
      subroutine Multiphase_thermalForcing(solnData, del, logc, higc)
        integer, dimension(MDIM),intent(IN) :: logc, higc 
        real, pointer, dimension(:, :, :, :) :: solnData
        real,dimension(MDIM), intent(IN) :: del
      end subroutine Multiphase_thermalForcing
   end interface

   interface
      subroutine Multiphase_divergence(solnData, facexData, faceyData, facezData, del, lo, hi)
        real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
        real,dimension(MDIM), intent(IN) :: del
        integer, dimension(MDIM),intent(IN) :: lo, hi
      end subroutine Multiphase_divergence
   end interface

   interface
      subroutine Multiphase_extrapFluxes(solnData,del,lo,hi, iteration)
        real, pointer, dimension(:, :, :, :) :: solnData
        real, dimension(MDIM), intent(IN) :: del
        integer, dimension(MDIM), intent(IN) :: lo,hi
        integer, intent(in) :: iteration
      end subroutine Multiphase_extrapFluxes
   end interface

   interface
      subroutine Multiphase_setMassFlux(solnData, del)
        real, pointer, dimension(:, :, :, :) :: solnData
        real, dimension(MDIM) :: del
      end subroutine Multiphase_setMassFlux
   end interface

   interface
      subroutine Multiphase_velForcing(solnData, facexData, faceyData, facezData, del, lo, hi, dt)
        implicit none
        real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
        real, dimension(MDIM), intent(in) :: del
        integer,dimension(MDIM), intent(in) :: lo,hi
        real, INTENT(IN) :: dt
      end subroutine Multiphase_velForcing
   end interface

   interface
      subroutine Multiphase_reInitGridVars(solnData)
        real, pointer, dimension(:, :, :, :) :: solnData
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
