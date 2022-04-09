!!***if* source/Simulation/SimulationForcing/incompFlow/Inlet/sim_inletVelBlk
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
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine sim_inletVelBlk2d(u, v, ru, rv, xcell, ycell, &
                             boundBox, dt, dx, dy, ix1, ix2, jy1, jy2, &
                             inletFlag, outletFlag, inletBuffer, inletGrowthRate, &
                             xMin, xMax, yMin, yMax, gravX, gravY)

   implicit none
   real, dimension(:, :, :), intent(in)    :: u, v
   real, dimension(:, :, :), intent(inout) :: ru, rv
   real, dimension(:), intent(in)          :: xcell, ycell
   real, dimension(:, :), intent(in)       :: boundBox
   real, intent(in)                        :: dt, dx, dy
   integer, intent(in)                     :: ix1, ix2, jy1, jy2
   integer, dimension(LOW:HIGH, MDIM), intent(in) :: inletFlag, outletFlag
   real, intent(in) :: inletBuffer, inletGrowthRate
   real, intent(in) :: xMin, xMax, yMin, yMax, gravX, gravY

end subroutine sim_inletVelBlk2d

subroutine sim_inletVelBlk3d(u, v, w, ru, rv, rw, xcell, ycell, zcell, &
                             boundBox, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, &
                             inletFlag, outletFlag, inletBuffer, inletGrowthRate, &
                             xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ)

   implicit none
   real, dimension(:, :, :), intent(in)  :: u, v, w
   real, dimension(:, :, :), intent(inout) :: ru, rv, rw
   real, dimension(:), intent(in)        :: xcell, ycell, zcell
   real, dimension(:, :), intent(in)     :: boundBox
   real, intent(in)                      :: dt, dx, dy, dz
   integer, intent(in)                   :: ix1, ix2, jy1, jy2, kz1, kz2
   integer, dimension(LOW:HIGH, MDIM), intent(in) :: inletFlag, outletFlag
   real, intent(in) :: inletBuffer, inletGrowthRate
   real, intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax, gravX, gravY, gravZ

end subroutine sim_inletVelBlk3d
