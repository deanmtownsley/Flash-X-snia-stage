!!***if* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterApplyBCToRegion
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

subroutine sim_heaterApplyBCToRegion(level, ivar, gridDataStruct, regionData, coordinates, regionSize, &
                                   guard, face, axis, secondDir, thirdDir)

   use Driver_interface, ONLY: Driver_abort
   use sim_heaterData, ONLY: sim_heaterType, sim_heaterInfo, sim_numHeaters
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

   if (face == HIGH) call Driver_abort('[sim_heaterApplyBCToRegion] not configured for face == HIGH')

   if (ivar == TEMP_VAR) then
      do k = 1, ke
         do j = 1, je
            do i = 1, guard
               do htr = 1, sim_numHeaters

                  heater => sim_heaterInfo(htr)

                  if (coordinates(i, j, k, IAXIS) .gt. heater%xMin .and. &
                      coordinates(i, j, k, IAXIS) .lt. heater%xMax .and. &
                      coordinates(i, j, k, KAXIS) .gt. heater%zMin .and. &
                      coordinates(i, j, k, KAXIS) .lt. heater%zMax) then

                     regionData(i, j, k, ivar) = 2*heater%wallTemp - regionData(offset - i, j, k, ivar)

                  end if

               end do
            end do
         end do
      end do

#ifdef MULTIPHASE_EVAPORATION
   else if (ivar == DFUN_VAR) then
      do k = 1, ke
         do j = 1, je
            do i = 1, guard
               do htr = 1, sim_numHeaters

                  heater => sim_heaterInfo(htr)

                  if (coordinates(i, j, k, IAXIS) .gt. heater%xMin .and. &
                      coordinates(i, j, k, IAXIS) .lt. heater%xMax .and. &
                      coordinates(i, j, k, KAXIS) .gt. heater%zMin .and. &
                      coordinates(i, j, k, KAXIS) .lt. heater%zMax) then

                     dynamicAngle = heater%rcdAngle

                     veli = regionData(guard + 1, j, k, VELX_VAR)*regionData(guard + 1, j, k, NRMX_VAR)
#if NDIM == MDIM
                     veli = veli + regionData(guard + 1, j, k, VELZ_VAR)*regionData(guard + 1, j, k, NRMZ_VAR)
#endif
                     if (veli .ge. 0.0) then
                        if (abs(veli) .le. heater%velContact) then
                           dynamicAngle = ((heater%advAngle - heater%rcdAngle)/(2*heater%velContact))*abs(veli) + &
                                          (heater%advAngle + heater%rcdAngle)/2.0d0
                        else
                           dynamicAngle = heater%advAngle
                        end if
                     end if

                     regionData(i, j, k, ivar) = regionData(offset - i, j, k, ivar) - &
                                                 del(axis)*cos(dynamicAngle*acos(-1.0)/180)
                  end if

               end do
            end do
         end do
      end do
#endif

   end if

end subroutine sim_heaterApplyBCToRegion
