!!****if* source/physics/IncompNS/IncompNSMain/vardens/IncompNS_corrector
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
!!***
!!REORDER(4): face[xyz]Data
!!REORDER(4): solnData

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

subroutine IncompNS_corrector(dt)

  use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getTileIterator, &
                               Grid_releaseTileIterator
  use Grid_tile,        ONLY : Grid_tile_t
  use Grid_iterator,    ONLY : Grid_iterator_t
  use ins_interface,    ONLY : ins_corrector_vardens
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_getNStep
  use IncompNS_data

  implicit none
  include "Flash_mpi.h"
  !-----Argument List-----
  real,    INTENT(IN) :: dt

!------------------------------------------------------------------------------------------
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
#if NDIM < MDIM
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData
  real, dimension(NFACE_VARS,1,1,1) :: facezData
#else
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
#endif
  integer TA(2),count_rate
  real*8  ET
  real del(MDIM)
  integer :: NStep
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor

!------------------------------------------------------------------------------------------
#if NDIM < MDIM
  nullify(solnData,facexData,faceyData)
#else
  nullify(solnData,facexData,faceyData,facezData)
#endif
 
  CALL SYSTEM_CLOCK(TA(1),count_rate)

  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())

     call itor%currentTile(tileDesc)
     call tileDesc%deltas(del)

     ! Get Index Limits:
     blkLimits   = tileDesc%limits
     blkLimitsGC = tileDesc%blkLimitsGC

     ! Point to blocks center and face vars:
     call tileDesc%getDataPtr(solnData,  CENTER)
     call tileDesc%getDataPtr(facexData, FACEX)
     call tileDesc%getDataPtr(faceyData, FACEY)

#if NDIM == 3
     call tileDesc%getDataPtr(facezData, FACEZ)
#endif

     ! update divergence-free velocities (not on block boundary)
     call ins_corrector_vardens(facexData(VELC_FACE_VAR,:,:,:),&
                           faceyData(VELC_FACE_VAR,:,:,:),&
                           facezData(VELC_FACE_VAR,:,:,:),&
                           !------------------------------!
                           facexData(SIGM_FACE_VAR,:,:,:),&
                           faceyData(SIGM_FACE_VAR,:,:,:),&
                           facezData(SIGM_FACE_VAR,:,:,:),&
                           !------------------------------!
                           facexData(PGN1_FACE_VAR,:,:,:),&
                           faceyData(PGN1_FACE_VAR,:,:,:),&
                           facezData(PGN1_FACE_VAR,:,:,:),&
                           !------------------------------!
                           facexData(PGN2_FACE_VAR,:,:,:),&
                           faceyData(PGN2_FACE_VAR,:,:,:),&
                           facezData(PGN2_FACE_VAR,:,:,:),&
                           !------------------------------!
                           facexData(RHOF_FACE_VAR,:,:,:),&
                           faceyData(RHOF_FACE_VAR,:,:,:),&
                           facezData(RHOF_FACE_VAR,:,:,:),&
                           !------------------------------!
                           solnData(PRES_VAR,:,:,:),&
                           ins_rhoGas,dt,del(DIR_X),del(DIR_Y),del(DIR_Z),&
                           GRID_ILO,GRID_IHI,&
                           GRID_JLO,GRID_JHI,&
                           GRID_KLO,GRID_KHI)

     ! Release pointers:
     call tileDesc%releaseDataPtr(solnData,  CENTER)
     call tileDesc%releaseDataPtr(facexData, FACEX)
     call tileDesc%releaseDataPtr(faceyData, FACEY)

#if NDIM ==3
     call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

     call itor%next()

  end do
  call Grid_releaseTileIterator(itor)  

  !------------------------------------------------------------------------------------------------------
  ! FILL GUARDCELLS FOR FINAL VELOCITIES AND PRESSURE:
  ! ---- ---------- --- ----- ---------- --- --------
  ! The pressure fill is used to compute distributed forces on
  ! immersed bodies.
  gcMask = .FALSE.
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! u
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! v
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! w
#endif
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

  CALL SYSTEM_CLOCK(TA(2),count_rate)
  ET=REAL(TA(2)-TA(1))/count_rate
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total IncompNS Corrector Time =',ET

  return
end subroutine IncompNS_corrector
