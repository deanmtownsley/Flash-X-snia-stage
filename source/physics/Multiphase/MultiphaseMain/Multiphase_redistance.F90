!!***if* source/physics/Multiphase/MultiphaseMain/Multiphase_redistance
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!!
!!
!!***
!!REORDER(4): solnData

#include "Simulation.h"
#include "constants.h"
#include "Multiphase.h"

subroutine Multiphase_redistance()

  use Multiphase_data
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use Driver_interface,    ONLY : Driver_getNStep
  use Grid_interface,      ONLY : Grid_fillGuardCells, Grid_getTileIterator, Grid_releaseTileIterator
  use Grid_tile,           ONLY : Grid_tile_t
  use Grid_iterator,       ONLY : Grid_iterator_t
  use Stencils_interface,  ONLY : Stencils_lsRedistance2d, Stencils_lsRedistance3d

!-----------------------------------------------------------------------------------------
  implicit none
  include "Flash_mpi.h"
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
  real, pointer, dimension(:,:,:,:) :: solnData
  integer :: ii,ierr,i,j,k
  real del(MDIM)
  real :: lsDT,minCellDiag
  real :: volSum,volSumAll
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor
  integer TA(2),count_rate
  real*8  ET

!-----------------------------------------------------------------------------------------
  nullify(solnData)

  CALL SYSTEM_CLOCK(TA(1),count_rate)

  !-----------------------------------------------------------------------!
  !-------STEP 1. COPY SOLUTION TO RDFN-----------------------------------!
  !-----------------------------------------------------------------------!  
  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData,  CENTER)
     solnData(RDFN_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)
     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)  

  !-----------------------------------------------------------------------!
  !-------STEP 2. LEVEL SET REDISTANCING----------------------------------!
  !-----------------------------------------------------------------------!   
  do ii = 1,mph_lsIt
     !------------------------------
     !- Level set redistancing 
     !------------------------------
     call Grid_getTileIterator(itor, nodetype=LEAF)
     do while(itor%isValid())
       call itor%currentTile(tileDesc)
       call tileDesc%getDataPtr(solnData,  CENTER)
       call tileDesc%deltas(del)
       minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.+del(DIR_Z)**2)
       lsDT = minCellDiag/5.0d0
#if NDIM < MDIM
        !--------------------------------------------
        ! Call DFUN re-initialization routine for 2D:
        !--------------------------------------------
        call Stencils_lsRedistance2d(solnData(DFUN_VAR,:,:,:), &
                                     solnData(RDFN_VAR,:,:,:), &
                                     lsDT,del(DIR_X),del(DIR_Y),  &
                                     GRID_ILO,GRID_IHI, &
                                     GRID_JLO,GRID_JHI)

#else
       call Stencils_lsRedistance3d(solnData(DFUN_VAR,:,:,:), &
                                    solnData(RDFN_VAR,:,:,:), &
                                    lsDT,del(DIR_X),del(DIR_Y),del(DIR_Z),&
                                    GRID_ILO,GRID_IHI, &
                                    GRID_JLO,GRID_JHI, &
                                    GRID_KLO,GRID_KHI)
#endif
        call tileDesc%releaseDataPtr(solnData, CENTER)
        call itor%next()
     end do
     call Grid_releaseTileIterator(itor)  

  !-----Fill distance function guard cells after each re-initialization to
  !communicate updates
  gcMask = .FALSE.
  gcMask(DFUN_VAR) = .TRUE.
  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
  !*********************************************************************************************************
  end do  ! End do: ii=1,lsit

  !-----------------------------------------------------------------------!
  !-------STEP 3. CALCULATE GAS VOLUME------------------------------------!
  !-----------------------------------------------------------------------!  
  volSum = 0.0
  volSumAll = 0.0
 
  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData,  CENTER)
     blkLimits = tileDesc%limits
     do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
      do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           if (solnData(RDFN_VAR,i,j,k) .gt. 0) then
#if NDIM < MDIM
           volSum = volSum + del(DIR_X) * del(DIR_Y)
#else
           volSum = volSum + del(DIR_X) * del(DIR_Y) * del(DIR_Z)
#endif
           end if
        end do
       end do
     end do
     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)  

  call MPI_Allreduce(volSum, volSumAll, 1, FLASH_REAL,&
                       MPI_SUM, MPI_COMM_WORLD, ierr)
  if (mph_meshMe .eq. 0) print*,"----------------------------------------"
  if (mph_meshMe .eq. 0) print*,"Total Gas Volume: ",volSumAll
  if (mph_meshMe .eq. 0) print*,"----------------------------------------"

  CALL SYSTEM_CLOCK(TA(2),count_rate)
  ET=REAL(TA(2)-TA(1))/count_rate
  if (mph_meshMe .eq. MASTER_PE)  write(*,*) 'Total Multiphase Redistance Time =',ET

  return

end subroutine Multiphase_redistance
