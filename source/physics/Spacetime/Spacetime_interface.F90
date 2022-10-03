!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief Public interfaces for Spacetime
!!
!! @details This is the header file for the Spacetime unit that defines
!!          its public interfaces.

!> @ingroup Spacetime
!! Interfaces to Spacetime public procedures
module Spacetime_interface

#include "constants.h"

   implicit none

   interface
      subroutine Spacetime_init
      end subroutine Spacetime_init
   end interface

   interface
      subroutine Spacetime_finalize
      end subroutine Spacetime_finalize
   end interface

   interface
      subroutine Spacetime(t, dt)
         implicit none
         real, intent(in) :: t, dt
      end subroutine Spacetime
   end interface

   interface
      subroutine Spacetime_computeDt(tileDesc, solnData, dtMin, dtMinLoc)
         use Grid_tile, only: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         real, pointer :: solnData(:, :, :, :)
         real, intent(inout) :: dtMin
         real, intent(inout) :: dtMinLoc(5)
      end subroutine Spacetime_computeDt
   end interface

   ! ADM accessors

   interface
      subroutine Spacetime_getLapse(alp, &
                                    tileDesc, solnData, loc)
         use Grid_tile, only: Grid_tile_t
         implicit none
         real, intent(out) :: alp
         type(Grid_tile_t), intent(in) :: tileDesc
         real, pointer :: solnData(:, :, :, :)
         integer, intent(in) :: loc(MDIM)
      end subroutine Spacetime_getLapse
   end interface

   interface
      subroutine Spacetime_getShift(betax, betay, betaz, &
                                    tileDesc, solnData, loc)
         use Grid_tile, only: Grid_tile_t
         implicit none
         real, intent(out) :: betax, betay, betaz
         type(Grid_tile_t), intent(in) :: tileDesc
         real, pointer :: solnData(:, :, :, :)
         integer, intent(in) :: loc(MDIM)
      end subroutine Spacetime_getShift
   end interface

   interface
      subroutine Spacetime_getMetric(gxx, gxy, gxz, gyy, gyz, gzz, &
                                     tileDesc, solnData, loc)
         use Grid_tile, only: Grid_tile_t
         implicit none
         real, intent(out) :: gxx, gxy, gxz, gyy, gyz, gzz
         type(Grid_tile_t), intent(in) :: tileDesc
         real, pointer :: solnData(:, :, :, :)
         integer, intent(in) :: loc(MDIM)
      end subroutine Spacetime_getMetric
   end interface

   interface
      subroutine Spacetime_getCurvature(Kxx, Kxy, Kxz, Kyy, Kyz, Kzz, &
                                        tileDesc, solnData, loc)
         use Grid_tile, only: Grid_tile_t
         implicit none
         real, intent(out) :: Kxx, Kxy, Kxz, Kyy, Kyz, Kzz
         type(Grid_tile_t), intent(in) :: tileDesc
         real, pointer :: solnData(:, :, :, :)
         integer, intent(in) :: loc(MDIM)
      end subroutine Spacetime_getCurvature
   end interface

   ! Stress-energy tensor updates

   interface
      subroutine Spacetime_addToStressEnergyTensor(E, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
                                                   tileDesc, solnData, loc)
         use Grid_tile, only: Grid_tile_t
         implicit none
         real, intent(in) :: E
         real, intent(in) :: Sx, Sy, Sz
         real, intent(in) :: Sxx, Sxy, Sxz, Syy, Syz, Szz
         type(Grid_tile_t), intent(in) :: tileDesc
         real, pointer :: solnData(:, :, :, :)
         integer, intent(in) :: loc(MDIM)
      end subroutine Spacetime_addToStressEnergyTensor
   end interface

   ! For constraint-based solvers

   interface
      subroutine Spacetime_applyHamiltonianConstraint(t)
         implicit none
         real, intent(in) :: t
      end subroutine Spacetime_applyHamiltonianConstraint
   end interface

   interface
      subroutine Spacetime_applyMomentumConstraint(t)
         implicit none
         real, intent(in) :: t
      end subroutine Spacetime_applyMomentumConstraint
   end interface

   ! For dynamic solvers

   interface Spacetime_molExplicitRHS
      subroutine Spacetime_molExplicitRHS(t, activeRHS)
         implicit none
         real, intent(in) :: t
         integer, intent(in) :: activeRHS
      end subroutine Spacetime_molExplicitRHS

      subroutine Spacetime_molExplicitRHS_tile(tileDesc, t, activeRHS)
         use Grid_tile, only: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         real, intent(in) :: t
         integer, intent(in) :: activeRHS
      end subroutine Spacetime_molExplicitRHS_tile
   end interface Spacetime_molExplicitRHS

   interface Spacetime_molFastRHS
      subroutine Spacetime_molFastRHS(t, activeRHS)
         implicit none
         real, intent(in) :: t
         integer, intent(in) :: activeRHS
      end subroutine Spacetime_molFastRHS

      subroutine Spacetime_molFastRHS_tile(tileDesc, t, activeRHS)
         use Grid_tile, only: Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(in) :: tileDesc
         real, intent(in) :: t
         integer, intent(in) :: activeRHS
      end subroutine Spacetime_molFastRHS_tile
   end interface Spacetime_molFastRHS

   interface
      subroutine Spacetime_molPreEvolve(t)
         implicit none
         real, intent(in) :: t
      end subroutine Spacetime_molPreEvolve
   end interface

   interface
      subroutine Spacetime_molPostUpdate(t)
         implicit none
         real, intent(in) :: t
      end subroutine Spacetime_molPostUpdate
   end interface

   interface
      subroutine Spacetime_molPostFastUpdate(t)
         implicit none
         real, intent(in) :: t
      end subroutine Spacetime_molPostFastUpdate
   end interface

   interface
      subroutine Spacetime_molPostTimeStep(t)
         implicit none
         real, intent(in) :: t
      end subroutine Spacetime_molPostTimeStep
   end interface

   interface
      subroutine Spacetime_molPostRegrid(t)
         implicit none
         real, intent(in) :: t
      end subroutine Spacetime_molPostRegrid
   end interface

   interface
      subroutine Spacetime_unitTest(fileUnit, perfect)
         implicit none
         integer, intent(in) :: fileUnit
         logical, intent(inout) :: perfect
      end subroutine Spacetime_unitTest
   end interface

end module Spacetime_interface
