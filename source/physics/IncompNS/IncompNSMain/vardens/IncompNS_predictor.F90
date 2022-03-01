!!****if* source/physics/IncompNS/IncompNSMain/vardens/IncompNS_predictor
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
!!***
!!REORDER(4): face[xyz]Data
!!REORDER(4): solnData

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

subroutine IncompNS_predictor(dt)

  use Grid_interface,      ONLY : Grid_fillGuardCells,Grid_getTileIterator,&
                                  Grid_releaseTileIterator
  use Grid_tile,           ONLY : Grid_tile_t
  use Grid_iterator,       ONLY : Grid_iterator_t
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use Driver_interface,    ONLY : Driver_getNStep
  use Stencils_interface,  ONLY : Stencils_integrateEuler
  use IncompNS_data

  implicit none
  include "Flashx_mpi.h"
  !-----Argument-List-----!
  real,    INTENT(IN) :: dt 

!------------------------------------------------------------------------------------------
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
#if NDIM < MDIM
  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData
  real, dimension(NFACE_VARS,1,1,1) :: facezData
#else
  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData
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
  !
  CALL SYSTEM_CLOCK(TA(1),count_rate)
  !------------------------------------------------------------------------------------------------------
  ! FILL GUARDCELLS:
  ! ---- ---------- --- ----- ---------- --- --------
  gcMask = .FALSE.
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! u
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! v
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! w
#endif
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
  !
  !------------------------------------------------------------------------------------------------------
  ! COMPUTE RIGHT HAND SIDE AND PREDICTOR STEP:
  ! ------- ----- ---- ---- --- --------- ----
  !------------------------------------------------------------------------------------------------------
  call Grid_getTileIterator(itor, nodetype=LEAF)
  !
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     blkLimits   = tileDesc%limits
     blkLimitsGC = tileDesc%blkLimitsGC
     call tileDesc%deltas(del)
     call tileDesc%getDataPtr(solnData,  CENTER)
     call tileDesc%getDataPtr(facexData, FACEX)
     call tileDesc%getDataPtr(faceyData, FACEY)

     call Stencils_integrateEuler(facexData(VELC_FACE_VAR,:,:,:),&
                                  facexData(RHDS_FACE_VAR,:,:,:),&
                                  dt,&
                                  GRID_ILO,GRID_IHI+1,&
                                  GRID_JLO,GRID_JHI,&
                                  GRID_KLO,GRID_KHI,&
                                  iSource=-ins_dpdx+ins_gravX)

     call Stencils_integrateEuler(faceyData(VELC_FACE_VAR,:,:,:),&
                                  faceyData(RHDS_FACE_VAR,:,:,:),&
                                  dt,&
                                  GRID_ILO,GRID_IHI,&
                                  GRID_JLO,GRID_JHI+1,&
                                  GRID_KLO,GRID_KHI,&
                                  iSource=-ins_dpdy+ins_gravY)

     ! Release pointers:
     call tileDesc%releaseDataPtr(solnData,  CENTER)
     call tileDesc%releaseDataPtr(facexData, FACEX)
     call tileDesc%releaseDataPtr(faceyData, FACEY)

#if NDIM ==3
     call tileDesc%getDataPtr(facezData, FACEZ)

     call Stencils_integrateEuler(facezData(VELC_FACE_VAR,:,:,:),&
                                  facezData(RHDS_FACE_VAR,:,:,:),&
                                  dt,&
                                  GRID_ILO,GRID_IHI,&
                                  GRID_JLO,GRID_JHI,&
                                  GRID_KLO,GRID_KHI+1,&
                                  iSource=-ins_dpdz+ins_gravZ)

     call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
  !------------------------------------------------------------------------------------------------------
  ! APPLY BC AND FILL GUARDCELLS FOR INTERMEDIATE VELOCITIES:
  ! ----- -- --- ---- ---------- --- ------------ ----------
  !------------------------------------------------------------------------------------------------------
  gcMask = .FALSE.
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! ustar
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! vstar
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! wstar
#endif
  ins_predcorrflg = .true.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
  ins_predcorrflg = .false.

  CALL SYSTEM_CLOCK(TA(2),count_rate)
  ET=REAL(TA(2)-TA(1))/count_rate
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total IncompNS Predictor Time =',ET

  return
end subroutine IncompNS_predictor
