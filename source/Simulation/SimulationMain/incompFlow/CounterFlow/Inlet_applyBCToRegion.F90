!!***if* source/Simulation/SimulationMain/incompFlow/CounterFlow/Inlet_applyBCToRegion
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
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine Inlet_applyBCToRegion(level, ivar, gridDataStruct, regionData, coordinates, regionSize, &
                                 guard, face, axis, secondDir, thirdDir)

   use Driver_interface, ONLY: Driver_getSimTime, Driver_abort

   use Simulation_data, ONLY: sim_channelDepth, sim_xMin, sim_xMax, &
                              sim_nozzleAmp, sim_nozzleFreq, sim_liqFlowRate, &
                              sim_gasFlowRate

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

   integer :: je, ke
   integer :: i, j, k, offset
   real, dimension(MDIM)  :: del
   real :: channelProfile, time
   logical :: isFace
   real, parameter :: pi = acos(-1.0)

   call Grid_getDeltas(level, del)
   call Driver_getSimTime(time)

   je = regionSize(SECOND_DIR)
   ke = regionSize(THIRD_DIR)

   isFace = (gridDataStruct == FACEX) .and. (axis == IAXIS)
   isFace = isFace .or. ((gridDataStruct == FACEY) .and. (axis == JAXIS))
   isFace = isFace .or. ((gridDataStruct == FACEZ) .and. (axis == KAXIS))

   if (face == LOW) then
      if (axis == JAXIS) then
         if (ivar == DFUN_VAR) then
            offset = 2*guard+1
            do k = 1, ke
               do j = 1, je
                  do i = 1, guard

                     !channelProfile = min(coordinates(i, j, k, IAXIS) - (sim_xMin + sim_channelDepth), &
                     !                     (sim_xMax - sim_channelDepth) - coordinates(i, j, k, IAXIS))

                     channelProfile = min(coordinates(i, j, k, IAXIS)-(sim_xMin+sim_channelDepth)- &
                                          sim_nozzleAmp*cos(sim_nozzleFreq*time*2*pi), &
                                          (sim_xMax-sim_channelDepth)-coordinates(i, j, k, IAXIS)+ &
                                          sim_nozzleAmp*cos(sim_nozzleFreq*time*2*pi+pi))

                     regionData(i, j, k, ivar) = 2*channelProfile-regionData(offset-i, j, k, ivar)

                  end do
               end do
            end do

         else if (ivar == VELC_FACE_VAR) then

            if (isFace) then
               offset = 2*guard+2
               do k = 1, ke
                  do j = 1, je
                     do i = 1, guard+1

                        !channelProfile = min(coordinates(i, j, k, IAXIS) - (sim_xMin + sim_channelDepth), &
                        !                     (sim_xMax - sim_channelDepth) - coordinates(i, j, k, IAXIS))

                        channelProfile = min(coordinates(i, j, k, IAXIS)-(sim_xMin+sim_channelDepth)- &
                                             sim_nozzleAmp*cos(sim_nozzleFreq*time*2*pi), &
                                             (sim_xMax-sim_channelDepth)-coordinates(i, j, k, IAXIS)+ &
                                             sim_nozzleAmp*cos(sim_nozzleFreq*time*2*pi+pi))

                        regionData(i, j, k, ivar) = ((1+sign(1., channelProfile))/2)*(sim_gasFlowRate)+ &
                                                    ((1-sign(1., channelProfile))/2)*(sim_liqFlowRate)

                     end do
                  end do
               end do
            end if

         end if
      end if
   else
      call Driver_abort('[Inlet_applyBCToRegion] not configured for face == HIGH')

   end if

end subroutine Inlet_applyBCToRegion
