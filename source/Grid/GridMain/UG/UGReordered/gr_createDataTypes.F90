!!****if* source/Grid/GridMain/UG/UGReordered/gr_createDataTypes
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
!!  gr_createDataTypes
!!
!! SYNOPSIS
!!
!!  call gr_createDataTypes()
!!  
!! DESCRIPTION 
!!  create the MPI derived Datatypes in the uniform grid for exchanging
!!  guardcells.  It is more efficient to make a datatype than to do multiple
!!  mpi send/recvs
!!  
!!  
!! ARGUMENTS 
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_createDataTypes()

  use Grid_data, ONLY : gr_exch, gr_gridDataStruct, gr_gridDataStructSize,&
       gr_numDataStruct
  use Grid_interface, ONLY :Grid_getBlkIndexLimits

implicit none
#include "constants.h"
#include "Simulation.h"
  include "Flashx_mpi.h"

  integer :: ierr,  size, stride, flashCont
  integer,dimension(MDIM)::blkExtent,guard
  integer :: exch1, exch2
  integer :: blockID=1
  integer,dimension(LOW:HIGH,MDIM):: blkLimits,blkLimitsGC

  integer :: i,j

  


  do i = 1,gr_numDataStruct
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,gr_gridDataStruct(i))
     blkExtent(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
     guard(:)=blkLimits(LOW,:)-blkLimitsGC(LOW,:)
  
     if(NDIM==1) then
        call MPI_TYPE_VECTOR(gr_gridDataStructSize(i),guard(IAXIS), blkExtent(IAXIS),&
             FLASH_REAL, gr_exch(i,IAXIS), ierr)  
        call MPI_TYPE_COMMIT(gr_exch(i,IAXIS), ierr)

        ! ------------------------ 2d -------------------------
     elseif(NDIM==2)then

        ! for up and down, ydir
        call MPI_TYPE_VECTOR(gr_gridDataStructSize(i), &
             blkExtent(IAXIS) * guard(JAXIS), &
             blkExtent(IAXIS) * blkExtent(JAXIS), &
             FLASH_REAL, gr_exch(i,JAXIS), ierr)
        call MPI_TYPE_COMMIT(gr_exch(i,JAXIS),ierr)
     
        ! for left and right, xdir
        call MPI_TYPE_VECTOR(blkExtent(JAXIS) * gr_gridDataStructSize(i), &
             guard(IAXIS), &
             blkExtent(IAXIS), FLASH_REAL, gr_exch(i,IAXIS), ierr)
        call MPI_TYPE_COMMIT(gr_exch(i,IAXIS), ierr)
        
     else
        ! ------------------------ 3d -------------------------
        ! for left and right, xdir
        call MPI_TYPE_VECTOR(blkExtent(JAXIS) * blkExtent(KAXIS) * &
             gr_gridDataStructSize(i), &
             guard(IAXIS), blkExtent(IAXIS), FLASH_REAL, gr_exch(i,IAXIS), ierr)
        call MPI_TYPE_COMMIT(gr_exch(i,IAXIS),ierr)

        !for up and down, ydir
        call MPI_TYPE_VECTOR(blkExtent(KAXIS) * gr_gridDataStructSize(i), &
             blkExtent(IAXIS) * guard(JAXIS), &
             blkExtent(IAXIS) * blkExtent(JAXIS), &
             FLASH_REAL, gr_exch(i,JAXIS), ierr)
        call MPI_TYPE_COMMIT(gr_exch(i,JAXIS),ierr)
        
        ! for back and forth, zdir
        call MPI_TYPE_VECTOR(gr_gridDataStructSize(i), &
             blkExtent(IAXIS) * blkExtent(JAXIS) * guard(KAXIS),  &
             blkExtent(IAXIS) * blkExtent(JAXIS) * blkExtent(KAXIS), &
             FLASH_REAL, gr_exch(i,KAXIS), ierr)
        call MPI_TYPE_COMMIT(gr_exch(i,KAXIS),ierr)
     end if
  end do
end subroutine gr_createDataTypes


