!!***if* source/physics/Multiphase/MultiphaseMain/Multiphase_setPressureJumps
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
!!REORDER(4): solnData, face[xyz]Data

#include "constants.h"
#include "Multiphase.h"
#include "Simulation.h"

subroutine Multiphase_setPressureJumps()

  use Multiphase_data
  use mph_interface,     ONLY : mph_setWeberJumps2d,mph_setWeberJumps3d
  use mph_evapInterface, ONLY : mph_setEvapJumps2d,mph_setEvapJumps3d
  use Timers_interface,  ONLY : Timers_start, Timers_stop
  use Driver_interface,  ONLY : Driver_getNStep
  use Grid_interface,    ONLY : Grid_fillGuardCells, Grid_getTileIterator, &
                                Grid_releaseTileIterator
  use Grid_tile,         ONLY : Grid_tile_t
  use Grid_iterator,     ONLY : Grid_iterator_t

!----------------------------------------------------------------------------------------------------------
  implicit none
  include "Flash_mpi.h"
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData
  integer :: ierr,i,j,k
  real del(MDIM)
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor
  integer TA(2),count_rate
  real*8  ET


!---------------------------------------------------------------------------------------------------------
  CALL SYSTEM_CLOCK(TA(1),count_rate)

  nullify(solnData,facexData,faceyData,facezData)

  !-------------------------------------------------
  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData,  CENTER)
     call tileDesc%getDataPtr(facexData,  FACEX)
     call tileDesc%getDataPtr(faceyData,  FACEY)
#if NDIM == MDIM
     call tileDesc%getDataPtr(facezData,  FACEZ)
#endif
     call tileDesc%deltas(del)

#if NDIM < MDIM
     call mph_setWeberJumps2d(solnData(DFUN_VAR,:,:,:),&
                              facexData(mph_iJumpVar,:,:,:),&
                              faceyData(mph_iJumpVar,:,:,:),&
                              del(DIR_X),del(DIR_Y),&
                              mph_invWeber,&
                              GRID_ILO_GC,GRID_IHI_GC, &
                              GRID_JLO_GC,GRID_JHI_GC)

#ifdef MULTIPHASE_EVAPORATION
     call mph_setEvapJumps2d(solnData(DFUN_VAR,:,:,:),&
                             facexData(mph_iJumpVar,:,:,:),&
                             faceyData(mph_iJumpVar,:,:,:),&
                             solnData(MFLX_VAR,:,:,:),mph_rhoGas,&
                             del(DIR_X),del(DIR_Y),&
                             GRID_ILO_GC,GRID_IHI_GC, &
                             GRID_JLO_GC,GRID_JHI_GC)
#endif

#else
     call mph_setWeberJumps3d(solnData(DFUN_VAR,:,:,:),&
                              facexData(mph_iJumpVar,:,:,:),&
                              faceyData(mph_iJumpVar,:,:,:),&
                              facezData(mph_iJumpVar,:,:,:),&
                              del(DIR_X),del(DIR_Y),del(DIR_Z),&
                              mph_invWeber,&
                              GRID_ILO_GC,GRID_IHI_GC, &
                              GRID_JLO_GC,GRID_JHI_GC, &
                              GRID_KLO_GC,GRID_KHI_GC)
#ifdef MULTIPHASE_EVAPORATION
     call mph_setEvapJumps3d(solnData(DFUN_VAR,:,:,:),&
                             facexData(mph_iJumpVar,:,:,:),&
                             faceyData(mph_iJumpVar,:,:,:),&
                             facezData(mph_iJumpVar,:,:,:),&
                             solnData(MFLX_VAR,:,:,:),mph_rhoGas,&
                             del(DIR_X),del(DIR_Y),del(DIR_Z),&
                             GRID_ILO_GC,GRID_IHI_GC, &
                             GRID_JLO_GC,GRID_JHI_GC, &
                             GRID_KLO_GC,GRID_KHI_GC)
#endif


#endif 
      ! Release pointers:
      call tileDesc%releaseDataPtr(solnData, CENTER)
      call tileDesc%releaseDataPtr(facexData, FACEX)
      call tileDesc%releaseDataPtr(faceyData, FACEY)
#if NDIM == MDIM
      call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif
      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)  

  gcMask = .FALSE.
  gcMask(NUNK_VARS+mph_iJumpVar) = .TRUE.
  gcMask(NUNK_VARS+1*NFACE_VARS+mph_iJumpVar) = .TRUE.
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+mph_iJumpVar) = .TRUE.
#endif
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

   CALL SYSTEM_CLOCK(TA(2),count_rate)
   ET=REAL(TA(2)-TA(1))/count_rate
   if (mph_meshMe .eq. MASTER_PE)  write(*,*) 'Multiphase setPressureJumps Time =',ET

   return

end subroutine Multiphase_setPressureJumps
