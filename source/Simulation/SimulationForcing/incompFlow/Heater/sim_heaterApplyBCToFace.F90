!!***if* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterApplyBCToFace
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
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine sim_heaterApplyBCToFace(level, ivar, gridDataStruct, regionData, coordinates, regionSize, &
                                   guard, face, axis, secondDir, thirdDir)

   use sim_heaterData
   use Grid_interface, ONLY: Grid_getDeltas

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

!-------------------------------------------------------------------------------------------
   integer :: je, ke
   integer :: i, j, k, htr, offset
   type(sim_heaterType), pointer :: heater
   real, dimension(MDIM)  :: del
   real :: dynamicAngle, veli

   call Grid_getDeltas(level, del)

   je = regionSize(SECOND_DIR)
   ke = regionSize(THIRD_DIR)

   offset = 2*guard + 1

end subroutine sim_heaterApplyBCToFace
