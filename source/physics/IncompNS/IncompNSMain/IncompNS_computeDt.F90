!!****if* source/physics/IncompNS/IncompNSMain/IncompNS_computeDt
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
!!  IncompNS_computeDt
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

subroutine IncompNS_computeDt(ins_mindt,ins_minloc)

  use Grid_interface,   ONLY : Grid_getTileIterator, &
                               Grid_releaseTileIterator, Grid_getNumBlksFromType, &
                               Grid_getCellVolumes

  use Grid_tile, ONLY : Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t
 
 
  use ins_interface, ONLY : ins_computeDtLocal

  use IncompNS_data, ONLY : ins_useIncompNS

  implicit none
#include "constants.h"
#include "Simulation.h"
  real, intent(INOUT) :: ins_mindt
  integer, intent(INOUT) :: ins_minloc(5)

  ! Local Variables:
  real, PARAMETER :: MAX_TSTEP = huge(1.0)
  real    :: dtLocal
  integer :: lminloc(5)

  !!prepatory data structures for passing coords to timestep routines
  real, dimension(MDIM) :: del


  !!arrays which hold the starting and ending indices of a block
  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC

  !!coordinate information to be passed into physics  
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  integer :: isize,jsize,ksize

  integer :: i, blockID

  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor

  if (.NOT. ins_useIncompNS) RETURN

  nullify(facexData,faceyData,facezData)

  !!Initialize all timestep variables.
  dtLocal    = MAX_TSTEP
  lminloc(:) = 0     
 
  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())

     call itor%currentTile(tileDesc)

     blkLimits   = tileDesc%limits
     blkLimitsGC = tileDesc%blkLimitsGC
     blockID     = tileDesc%level

     call tileDesc%deltas(del)

     call tileDesc%getDataPtr(facexData, FACEX)
     call tileDesc%getDataPtr(faceyData, FACEY)

#if NDIM == 3
     call tileDesc%getDataPtr(facezData, FACEZ)
#endif

     isize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jsize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     ksize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifdef DEBUG_DRIVER
     print*,'going to call INS timestep'
#endif

     ! ins_computeDtLocal computed de blockID min_dt
     ! which if is .lt. dtLocal, is assigned to dtLocal
     call ins_computeDtLocal (blockID,         &
                         isize, jsize, ksize,  &
              del(IAXIS),del(JAXIS),del(KAXIS),&
                         blkLimits,blkLimitsGC,&
                         facexData,faceyData,  &
                         facezData,            &
                         dtLocal,lminloc)


     call tileDesc%releaseDataPtr(facexData, FACEX)
     call tileDesc%releaseDataPtr(faceyData, FACEY)

#if NDIM ==3
     call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

     call itor%next()

#ifdef DEBUG_DRIVER
     print*,'returned from INS timestep'
#endif

  end do
  call Grid_releaseTileIterator(itor)  


  ! Assign to  ins_mindt
  ins_mindt = dtLocal
  ins_minloc= lminloc
  return

end subroutine IncompNS_computeDt
