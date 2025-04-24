!!****f* source/physics/sourceTerms/Heater/HeaterMain/Heater_checkSites
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
!!***
!!REORDER(4): solnData, face[xyz]Data

#include "Simulation.h"
#include "constants.h"

subroutine Heater_checkSites(solnData,del,xcell,ycell,zcell,bndBox,logc,higc, lblock)

  use Timers_interface, ONLY: Timers_start, Timers_stop
  
  use Heater_data
  use Heater_type, ONLY: Heater_type_t
  use Driver_interface, ONLY: Driver_abort
  
  implicit none
  
  real, pointer, dimension(:, :, :, :) :: solnData
  real, dimension(MDIM),intent(IN) :: del
  real, dimension(:), intent(in)          :: xcell, ycell, zcell
  real, dimension(LOW:HIGH,MDIM), intent(in)        :: bndBox
  integer, dimension(MDIM),intent(in) :: logc, higc
  integer, intent(in) :: lblock
  
  
  type(Heater_type_t), pointer :: heater
  integer :: i, j, k, isite, htr, isiteblk
  real    :: xi, xp, yi, yp, zi, zp
  real    :: phiFSW, phiFSE, phiFNW, phiFNE
  real    :: phiBSW, phiBSE, phiBNW, phiBNE
  real    :: phiSW, phiSE,phiNW, phiNE
              
  real    :: phiSite
  integer :: ix1, ix2, jy1, jy2, kz1, kz2  
  
  !----------------------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------------------
  
  call Timers_start("Heater_checkSites")
  ix1=logc(IAXIS); ix2=higc(IAXIS)
  jy1=logc(JAXIS); jy2=higc(JAXIS)
  kz1=logc(KAXIS); kz2=higc(KAXIS)    
  
  
  
  do htr = 1, htr_numHeaters
     
     heater => htr_heaterInfo(htr)
     
     if (bndBox(HIGH, IAXIS) .le. heater%xMin .or. bndBox(LOW, IAXIS) .ge. heater%xMax .or. &
          bndBox(HIGH, JAXIS) .le. heater%yMin .or. bndBox(LOW, JAXIS) .ge. heater%yMax) cycle
#if NDIM < MDIM
     k = 1
     do j = jy1, jy2-1
        do i = ix1, ix2-1
           do isiteblk = 1, heater%numSitesBlk(lblock)
              isite = heater%siteMapOnProc(lblock, isiteblk)
              
              if (isite < 1) call Driver_abort("[htr_checkSitesBlk2d] isite < 1")
              
              xi = xcell(i)
              xp = xcell(i+1)
              yi = ycell(j)
              yp = ycell(j+1)
              
              phiSW = solnData(DFUN_VAR,i, j, k)
              phiSE = solnData(DFUN_VAR,i+1, j, k)
              phiNW = solnData(DFUN_VAR,i, j+1, k)
              phiNE = solnData(DFUN_VAR,i+1, j+1, k)
              
              phiSite = (phiSW+phiSE+phiNW+phiNE)/4.
              
              if (xi .le. heater%xSiteProc(isite) .and. xp .ge. heater%xSiteProc(isite) .and. &
                   yi .le. heater%ySiteProc(isite) .and. yp .ge. heater%ySiteProc(isite)) then
                 
                 if (phiSite .ge. 0.0) then
                    heater%siteIsAttachedCurr(isite) = heater%siteIsAttachedCurr(isite) .or. .true.
                 else
                    heater%siteIsAttachedCurr(isite) = heater%siteIsAttachedCurr(isite) .or. .false.
                 end if
              end if
           end do
        end do
     end do
#else
     if (bndBox(HIGH, IAXIS) .le. heater%xMin .or. bndBox(LOW, IAXIS) .ge. heater%xMax .or. &
          bndBox(HIGH, JAXIS) .le. heater%yMin .or. bndBox(LOW, JAXIS) .ge. heater%yMax .or. &
          bndBox(HIGH, KAXIS) .le. heater%zMin .or. bndBox(LOW, KAXIS) .ge. heater%zMax) cycle
     
     do k = kz1, kz2-1
        do j = jy1, jy2-1
           do i = ix1, ix2-1
              do isiteblk = 1, heater%numSitesBlk(lblock)
                 isite = heater%siteMapOnProc(lblock, isiteblk)
                 
                 if (isite < 1) call Driver_abort("[htr_checkSitesBlk3d] isite < 1")
                 
                 xi = xcell(i)
                 xp = xcell(i+1)
                 yi = ycell(j)
                 yp = ycell(j+1)
                 zi = zcell(k)
                 zp = zcell(k+1)
                 
                 phiFSW = solnData(DFUN_VAR,i, j, k)
                 phiFSE = solnData(DFUN_VAR,i+1, j, k)
                 phiFNW = solnData(DFUN_VAR,i, j+1, k)
                 phiFNE = solnData(DFUN_VAR,i+1, j+1, k)
                 
                 phiBSW = solnData(DFUN_VAR,i, j, k+1)
                 phiBSE = solnData(DFUN_VAR,i+1, j, k+1)
                 phiBNW = solnData(DFUN_VAR,i, j+1, k+1)
                 phiBNE = solnData(DFUN_VAR,i+1, j+1, k+1)
                 
                 phiSite = (phiFSW+phiFSE+phiFNW+phiFNE+phiBSW+phiBSE+phiBNW+phiBNE)/8.
                 
                 if (xi .le. heater%xSiteProc(isite) .and. xp .ge. heater%xSiteProc(isite) .and. &
                      yi .le. heater%ySiteProc(isite) .and. yp .ge. heater%ySiteProc(isite) .and. &
                      zi .le. heater%zSiteProc(isite) .and. zp .ge. heater%zSiteProc(isite)) then
                    
                    if (phiSite .ge. 0.0) then
                       heater%siteIsAttachedCurr(isite) = heater%siteIsAttachedCurr(isite) .or. .true.
                    else
                       heater%siteIsAttachedCurr(isite) = heater%siteIsAttachedCurr(isite) .or. .false.
                    end if
                 end if
              end do
           end do
        end do
     end do
#endif
  end do
  call Timers_stop("Heater_checkSites")
  
  return
end subroutine Heater_checkSites
