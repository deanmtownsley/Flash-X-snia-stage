!!****if* source/physics/IncompNS/IncompNSMain/constdens/ins_computeDtLocal
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
!!  ins_computeDtLocal
!!
!!
!! SYNOPSIS
!!
!!  
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!***

!! REORDER(4): facexData, faceyData, facezData

subroutine ins_computeDtLocal(blockID,   & 
                              isize, jsize, ksize,  &
                              dx, dy, dz,           &
                              blkLimits,blkLimitsGC,&
                              facexData,faceyData,  &
                              facezData,            &
                              dtLocal, lminloc )

  use IncompNS_data, ONLY : ins_cflflg, ins_cfl, ins_sigma, ins_invReynolds, ins_dtspec
  use IncompNS_data, ONLY : ins_meshMe

  implicit none

#include "Simulation.h"
#include "constants.h"
  integer, intent(IN) :: blockID
  integer,dimension(2,MDIM), intent(IN) :: blkLimits,blkLimitsGC
  integer, intent(IN) :: isize,jsize,ksize
  real, intent(IN) :: dx, dy, dz
  real, pointer,dimension(:,:,:,:)  :: facexData,faceyData,facezData
  real, intent(INOUT) :: dtLocal
  integer, intent(INOUT) :: lminloc(5)

  ! Local variables:
  real, parameter :: eps = 1.e-12
  real :: dtc,dtv,dtl,velcoeff


  if (ins_cflflg .eq. 0) then
     dtlocal    = ins_dtspec
     lminloc(:) = 0
     return
  endif

# if NDIM == MDIM

  velcoeff =  MAX( MAXVAL(ABS(facexData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
                                                      blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                                                      blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))/dx), &

                   MAXVAL(ABS(faceyData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                                                      blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,&
                                                      blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))/dx), &

                   MAXVAL(ABS(facezData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                                                      blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                                                      blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1))/dx))

  if (velcoeff .gt. eps) then
  dtc = ins_cfl / velcoeff
  else
  dtc = ins_cfl / eps
  endif
  
  dtv = ins_sigma / (ins_invReynolds*MAX( 2./(dx*dx), 2./(dy*dy),2./(dz*dz) ))
         
# elif NDIM == 2

  velcoeff =  MAX( MAXVAL(ABS(facexData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
                                                      blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),:))/dx), &
                   MAXVAL(ABS(faceyData(VELC_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                                                      blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1,:))/dx))

  if (velcoeff .gt. eps) then
  dtc = ins_cfl / velcoeff
  else
  dtc = ins_cfl / eps  
  endif
  
  dtv = ins_sigma / (ins_invReynolds*MAX( 1.0/(dx*dx), 1.0/(dy*dy)))

# endif

  dtl = MIN(dtc,dtv)

  if (dtl .lt. dtLocal) then
     dtLocal = dtl
     ! Cell located at center of Block - Used to define block location. 
     lminloc(IAXIS) = NGUARD + NXB/2
     lminloc(JAXIS) = NGUARD + NYB/2
#if NDIM == MDIM
     lminloc(KAXIS) = NGUARD + NZB/2
#else
     lminloc(KAXIS) = CONSTANT_ONE
#endif
     lminloc(4) = blockID
     lminloc(5) = ins_meshMe
  endif

  return

end subroutine ins_computeDtLocal
