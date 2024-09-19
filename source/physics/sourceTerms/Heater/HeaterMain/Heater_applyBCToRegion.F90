!!***if* source/physics/sourceTerms/Heater/HeaterMain/Heater_applyBCToRegion
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

subroutine Heater_applyBCToRegion(level, ivar, gridDataStruct, regionData, coordinates, regionSize, &
                                  guard, face, axis, secondDir, thirdDir)

   use Simulation_data, ONLY: sim_yMin, sim_yMax
   use Multiphase_data, ONLY: mph_thcoGas   
   use Driver_interface, ONLY: Driver_abort
   use Heater_data, ONLY: htr_heaterInfo, htr_numHeaters
   use Heater_type, ONLY: Heater_type_t
   use Grid_interface, ONLY: Grid_getDeltas
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get

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
   type(Heater_type_t), pointer :: heater
   real, dimension(MDIM)  :: del
   real :: dynamicAngle, veli
   real :: wall_dist_fun, th_con
   call Grid_getDeltas(level, del)

   call RuntimeParameters_get('ymin', sim_yMin)
   call RuntimeParameters_get('ymax', sim_yMax)
   call RuntimeParameters_get("mph_thcoGas", mph_thcoGas)

   th_con = mph_thcoGas

   je = regionSize(SECOND_DIR)
   ke = regionSize(THIRD_DIR)

   offset = 2*guard+1

!   if (face == HIGH) call Driver_abort('[Heater_applyBCToRegion] not configured for face == HIGH')

   if (face == HIGH) then
           if (ivar == TEMP_VAR) then
              do k = 1, ke
                 do j = 1, je
                    do i = 1, guard
                       do htr = 1, htr_numHeaters

                          heater => htr_heaterInfo(htr)

                          if (coordinates(i, j, k, IAXIS) .gt. heater%xMin .and. &
                              coordinates(i, j, k, IAXIS) .lt. heater%xMax .and. &
                              coordinates(i, j, k, KAXIS) .gt. heater%zMin .and. &
                              coordinates(i, j, k, KAXIS) .lt. heater%zMax) then

                             if(abs(0.5*(heater%yMax + heater%yMin) - sim_yMin) .gt. &
                                abs(0.5*(heater%yMax + heater%yMin) - sim_yMax))then

                                     if(heater%heat_flux_flag .eq. 1)then

                                          wall_dist_fun  = 0.5*(regionData(guard, j, k, DFUN_VAR)+&
                                                                 regionData(guard+1, j, k, DFUN_VAR))
                                          
                                          if(wall_dist_fun .le. 0) th_con = 1.0

                                          regionData(offset - i, j, k, ivar) = &
                                                                  ((heater%nd_heat_flux)*del(axis))*(offset - 2*i)/(th_con) +&
                                                                  regionData(i, j, k, ivar)

                                          !regionData(guard+i, j, k, ivar) = &
                                          !                        ((heater%nd_heat_flux)*del(axis))/(th_con) +&
                                          !                        regionData(guard+i-1, j, k, ivar)
                                          !regionData(guard+i, j, k, ivar) = &
                                          !                        0.33*(((heater%nd_heat_flux)*4*del(axis))/(th_con) +&
                                          !                        regionData(guard+i-1, j, k, ivar) +&
                                          !                        3*regionData(guard+i-2, j, k, ivar)-&
                                          !                        regionData(guard+i-3, j, k, ivar))
                                                                

                                     
                                     else
                                          if(heater%non_uniform_temp_flag .eq. 1)then

                                                regionData(offset - i, j, k, ivar) = &
                                                               2*(heater%C3*(coordinates(i, j, k, IAXIS))**3+&
                                                               heater%C2*(coordinates(i, j, k, IAXIS))**2+&
                                                               heater%C1*(coordinates(i, j, k, IAXIS))+&
                                                               heater%C0)-&
                                                               regionData(i, j, k, ivar)
                                          else        

                                                regionData(offset - i, j, k, ivar) = 2*heater%wallTemp-&
                                                                  regionData(i, j, k, ivar)    
                                          end if

                                     end if
                             end if

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
                       do htr = 1, htr_numHeaters

                          heater => htr_heaterInfo(htr)

                          if (coordinates(i, j, k, IAXIS) .gt. heater%xMin .and. &
                              coordinates(i, j, k, IAXIS) .lt. heater%xMax .and. &
                              coordinates(i, j, k, KAXIS) .gt. heater%zMin .and. &
                              coordinates(i, j, k, KAXIS) .lt. heater%zMax) then

                             dynamicAngle = heater%rcdAngle

                             veli = regionData(guard, j, k, VELX_VAR)*regionData(guard, j, k, NRMX_VAR)
#if NDIM == MDIM
                             veli = veli+regionData(guard, j, k, VELZ_VAR)*regionData(guard, j, k, NRMZ_VAR)
#endif
                             if (veli .ge. 0.0) then
                                if (abs(veli) .le. heater%velContact) then
                                   dynamicAngle = ((heater%advAngle-heater%rcdAngle)/(2*heater%velContact))*abs(veli)+ &
                                                  (heater%advAngle+heater%rcdAngle)/2.0d0
                                else
                                   dynamicAngle = heater%advAngle
                                end if
                             end if

                             if(abs(0.5*(heater%yMax + heater%yMin) - sim_yMin).gt. &
                                abs(0.5*(heater%yMax + heater%yMin) - sim_yMax))then
                             
                                     regionData(offset - i, j, k, ivar) = regionData(i, j, k, ivar)- &
                                                         del(axis)*cos(dynamicAngle*acos(-1.0)/180)
                                   
                                     !regionData(offset - i, j, k, ivar) = regionData(i, j, k, ivar)- &
                                     !                    del(axis)*(offset - 2*i)*cos(dynamicAngle*acos(-1.0)/180)                                                 
                                                 
                             end if

                          end if

                       end do
                    end do
                 end do
              end do
#endif

           end if

   else  
           if (ivar == TEMP_VAR) then
              do k = 1, ke
                 do j = 1, je
                    do i = 1, guard
                       do htr = 1, htr_numHeaters

                          heater => htr_heaterInfo(htr)

                          if (coordinates(i, j, k, IAXIS) .gt. heater%xMin .and. &
                              coordinates(i, j, k, IAXIS) .lt. heater%xMax .and. &
                              coordinates(i, j, k, KAXIS) .gt. heater%zMin .and. &
                              coordinates(i, j, k, KAXIS) .lt. heater%zMax) then

                             if(abs(0.5*(heater%yMax + heater%yMin) - sim_yMin) .lt. & 
                                abs(0.5*(heater%yMax + heater%yMin) - sim_yMax))then

                                     if(heater%heat_flux_flag .eq. 1)then

                                          wall_dist_fun  = 0.5*(regionData(guard, j, k, DFUN_VAR) +&
                                                           regionData(guard+1, j, k, DFUN_VAR))
                                          !write(*,*) regionData(guard+1,j,k,DFUN_VAR)
                                          !write(*,*) wall_dist_fun         

                                          if(wall_dist_fun .le. 0) th_con = 1.0

                                          !write(*,*) heater%nd_heat_flux
                                          !write(*,*) del(axis)
                                          !write(*,*) th_con

                                          regionData(i, j, k, ivar) = &
                                                                  ((heater%nd_heat_flux)*del(axis))*(offset-2*i)/(th_con) +&
                                                                  regionData(offset - i, j, k, ivar)

                                          !regionData(guard+1-i, j, k, ivar) = &
                                          !                        ((heater%nd_heat_flux)*del(axis))/(th_con) +&
                                          !                        regionData(guard+2 - i, j, k, ivar)

                                          !regionData(guard+1-i, j, k, ivar) = &
                                          !                        0.33*(((heater%nd_heat_flux)*4*del(axis))/(th_con) +&
                                          !                        regionData(guard+2-i, j, k, ivar) +&
                                          !                        3*regionData(guard+3-i, j, k, ivar)-&
                                          !                        regionData(guard+4-i, j, k, ivar))


                                     else  


                                          if(heater%non_uniform_temp_flag .eq. 1)then

                                                 regionData(i, j, k, ivar) = &
                                                               2*(heater%C3*(coordinates(i, j, k, IAXIS))**3+&
                                                               heater%C2*(coordinates(i, j, k, IAXIS))**2+&
                                                               heater%C1*(coordinates(i, j, k, IAXIS))+&
                                                               heater%C0)-&
                                                               regionData(offset - i, j, k, ivar)
                                          else

                                                 regionData(i, j, k, ivar) = 2*heater%wallTemp-&
                                                                  regionData(offset - i, j, k, ivar)
                                          end if

                                     end if
                             end if

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
                       do htr = 1, htr_numHeaters

                          heater => htr_heaterInfo(htr)

                          if (coordinates(i, j, k, IAXIS) .gt. heater%xMin .and. &
                              coordinates(i, j, k, IAXIS) .lt. heater%xMax .and. &
                              coordinates(i, j, k, KAXIS) .gt. heater%zMin .and. &
                              coordinates(i, j, k, KAXIS) .lt. heater%zMax) then

                             dynamicAngle = heater%rcdAngle

                             veli = regionData(guard+1, j, k, VELX_VAR)*regionData(guard+1, j, k, NRMX_VAR)
#if NDIM == MDIM
                             veli = veli+regionData(guard+1, j, k, VELZ_VAR)*regionData(guard+1, j, k, NRMZ_VAR)
#endif
                             if (veli .ge. 0.0) then
                                if (abs(veli) .le. heater%velContact) then
                                   dynamicAngle = ((heater%advAngle-heater%rcdAngle)/(2*heater%velContact))*abs(veli)+ &
                                                  (heater%advAngle+heater%rcdAngle)/2.0d0
                                else
                                   dynamicAngle = heater%advAngle
                                end if
                             end if

                             if(abs(0.5*(heater%yMax + heater%yMin) - sim_yMin) .lt. &
                                abs(0.5*(heater%yMax + heater%yMin) - sim_yMax))then

                                     regionData(i, j, k, ivar) = regionData(offset-i, j, k, ivar)- &
                                                           del(axis)*cos(dynamicAngle*acos(-1.0)/180)
                                  !   regionData(i, j, k, ivar) = regionData(offset-i, j, k, ivar)- &
                                  !                         del(axis)*(offset-2*i)*cos(dynamicAngle*acos(-1.0)/180)                                                   
    
                             end if
                          end if

                       end do
                    end do
                 end do
              end do
#endif

           end if

   end if        
end subroutine Heater_applyBCToRegion
