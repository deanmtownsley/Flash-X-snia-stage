!!****f* source/Simulation/SimulationMain/incompFlow/PoolBoiling/localAPI/heater/sim_heaterApplyBC
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

subroutine sim_heaterApplyBC(dt)

  use Simulation_data,     ONLY : sim_meshMe
  use Grid_interface,      ONLY : Grid_getTileIterator,Grid_releaseTileIterator,&
                                  Grid_fillGuardCells,Grid_getCellCoords
  use Grid_tile,           ONLY : Grid_tile_t
  use Grid_iterator,       ONLY : Grid_iterator_t
  use sim_heaterInterface, ONLY : sim_heaterApplyBCToBlk2d

  implicit none
  real, intent(in) :: dt

!----------------------------------------------------------------------------------------
  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
  integer, dimension(2,MDIM)        :: blkLimits, blkLimitsGC
  integer, dimension(MDIM)          :: lo, hi
  real, dimension(GRID_IHI_GC)      :: xCenter
  real, dimension(GRID_JHI_GC)      :: yCenter
  real, dimension(GRID_KHI_GC)      :: zCenter
  integer :: TA(2),count_rate
  real*8  :: ET
  real    :: del(MDIM)
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor
  real    :: boundBox(LOW:HIGH, 1:MDIM)

!----------------------------------------------------------------------------------------
  nullify(solnData,facexData,faceyData,facezData)

  CALL SYSTEM_CLOCK(TA(1),count_rate)

  call Grid_getTileIterator(itor, nodetype=LEAF)
  !
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     blkLimits   = tileDesc%limits
     blkLimitsGC = tileDesc%blkLimitsGC
     call tileDesc%deltas(del)
     call tileDesc%boundBox(boundBox)
     call tileDesc%getDataPtr(solnData,  CENTER)
     call tileDesc%getDataPtr(facexData, FACEX)
     call tileDesc%getDataPtr(faceyData, FACEY)

     lo=blkLimitsGC(LOW,:)
     hi=blkLimitsGC(HIGH,:)

     xCenter = 0.0
     yCenter = 0.0
     zCenter = 0.0
     call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xCenter)
     call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, lo, hi, yCenter)
     if (NDIM == MDIM) call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, lo, hi, zCenter)

     call sim_heaterApplyBCToBlk2d(solnData(DFRC_VAR,:,:,:),&
                                   solnData(TFRC_VAR,:,:,:),&
                                   solnData(DFUN_VAR,:,:,:),&
                                   solnData(TEMP_VAR,:,:,:),&
                                   xCenter,yCenter,boundBox,&
                                   dt,del(IAXIS),del(JAXIS),&
                                   GRID_ILO_GC,GRID_IHI_GC,&
                                   GRID_JLO_GC,GRID_JHI_GC)

#if NDIM == MDIM
     call tileDesc%getDataPtr(facezData, FACEZ)

     call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

     ! Release pointers:
     call tileDesc%releaseDataPtr(solnData,  CENTER)
     call tileDesc%releaseDataPtr(facexData, FACEX)
     call tileDesc%releaseDataPtr(faceyData, FACEY)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

   gcMask = .false.
   gcMask(DFUN_VAR) = .true.
   gcMask(TEMP_VAR) = .true.
   call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
        maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

  CALL SYSTEM_CLOCK(TA(2),count_rate)
  ET=REAL(TA(2)-TA(1))/count_rate
  if (sim_meshMe .eq. MASTER_PE)  write(*,*) 'Total sim_heater ApplyBC Time =',ET

  return

end subroutine sim_heaterApplyBC
