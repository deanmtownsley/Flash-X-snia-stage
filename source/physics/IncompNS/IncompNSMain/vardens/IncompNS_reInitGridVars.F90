!!****if* source/physics/IncompNS/IncompNSMain/vardens/IncompNS_reInitGridVars
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
!!
!! IncompNS_reInitGridVars
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!! 
!!REORDER(4): face[xyz]Data
!!REORDER(4): solnData

#include "constants.h"
#include "IncompNS.h"
#include "Simulation.h"

subroutine IncompNS_reInitGridVars()

  use Grid_interface,   ONLY : Grid_getTileIterator, Grid_releaseTileIterator
  use Grid_tile,        ONLY : Grid_tile_t
  use Grid_iterator,    ONLY : Grid_iterator_t
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_getNStep
  use IncompNS_data

!------------------------------------------------------------------------------------------
  implicit none
  include "Flashx_mpi.h"
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  integer TA(2),count_rate
  real*8  ET
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor

!------------------------------------------------------------------------------------------
  nullify(solnData,facexData,faceyData,facezData)

  CALL SYSTEM_CLOCK(TA(1),count_rate)

  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData,  CENTER)
     call tileDesc%getDataPtr(facexData, FACEX)
     call tileDesc%getDataPtr(faceyData, FACEY)
     solnData(VISC_VAR,:,:,:) = 1.
     facexData(RHOF_FACE_VAR,:,:,:) = 1.
     faceyData(RHOF_FACE_VAR,:,:,:) = 1.
     facexData(SIGM_FACE_VAR,:,:,:) = 0.
     faceyData(SIGM_FACE_VAR,:,:,:) = 0.
     facexData(RHDS_FACE_VAR,:,:,:) = 0.
     faceyData(RHDS_FACE_VAR,:,:,:) = 0.
     ! Release pointers:
     call tileDesc%releaseDataPtr(solnData,  CENTER)
     call tileDesc%releaseDataPtr(facexData, FACEX)
     call tileDesc%releaseDataPtr(faceyData, FACEY)
#if NDIM ==3
     call tileDesc%getDataPtr(facezData,  FACEZ)
     facezData(RHOF_FACE_VAR,:,:,:) = 1.
     facezData(SIGM_FACE_VAR,:,:,:) = 0.
     facezData(RHDS_FACE_VAR,:,:,:) = 0.
     call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)  

  CALL SYSTEM_CLOCK(TA(2),count_rate)
  ET=REAL(TA(2)-TA(1))/count_rate
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total IncompNS reInit Grid Vars Time =',ET

  return
end subroutine IncompNS_reInitGridVars
