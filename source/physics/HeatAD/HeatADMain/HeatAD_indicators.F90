!!***if* source/physics/HeatAD/HeatADMain/HeatAD_indicators
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
!!
!!***
!!REORDER(4): solnData

#include "constants.h"
#include "HeatAD.h"
#include "Simulation.h"   

subroutine HeatAD_indicators()

   use HeatAD_data
   use ht_interface,        ONLY : ht_indicators
   use Timers_interface,    ONLY : Timers_start, Timers_stop
   use Driver_interface,    ONLY : Driver_getNStep
   use Grid_interface,      ONLY : Grid_getTileIterator,Grid_releaseTileIterator
   use Grid_tile,           ONLY : Grid_tile_t
   use Grid_iterator,       ONLY : Grid_iterator_t

!--------------------------------------------------------------------------------------------
   implicit none
   include"Flashx_mpi.h"
   real ::  del(MDIM)
   integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:,:,:,:) :: solnData
   integer TA(2),count_rate,ierr
   real*8  ET
   type(Grid_tile_t) :: tileDesc
   type(Grid_iterator_t) :: itor
   real :: tempmaxaux, tempminaux, tempmax, tempmin

!---------------------------------------------------------------------------------------------
  nullify(solnData)

  !------------------------------------------------------------!
  !-------COMPUTE INDICATORS-----------------------------------!
  !------------------------------------------------------------!   
  tempmaxaux = -10.**(10.)
  tempminaux =  10.**(10.)
  call Grid_getTileIterator(itor, nodetype=LEAF)
   do while(itor%isValid())
       call itor%currentTile(tileDesc)
       call tileDesc%getDataPtr(solnData,  CENTER)

       call ht_indicators(solnData(TEMP_VAR,:,:,:),&
                          GRID_ILO,GRID_IHI,&
                          GRID_JLO,GRID_JHI,&
                          GRID_KLO,GRID_KHI,&
                          tempminaux,tempmaxaux)

       call tileDesc%releaseDataPtr(solnData, CENTER)
       call itor%next()
  end do
  call Grid_releaseTileIterator(itor)  

  call MPI_Allreduce(tempmaxaux, tempmax, 1, FLASH_REAL,&
                     MPI_MAX, ht_meshComm, ierr)

  call MPI_Allreduce(tempminaux, tempmin, 1, FLASH_REAL,&
                     MPI_MIN, ht_meshComm, ierr)

  if (ht_meshMe .eq. MASTER_PE) then
     write(*,*) ' '
     write(*,'(A24,2g14.6)') ' Min , Max  T =',tempmin,tempmax
  endif

end subroutine HeatAD_indicators
