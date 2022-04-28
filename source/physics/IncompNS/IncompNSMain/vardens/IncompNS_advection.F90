!!****if* source/physics/IncompNS/IncompNSMain/vardens/IncompNS_advection
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

subroutine IncompNS_advection()

  use Grid_interface,      ONLY : Grid_getTileIterator,Grid_releaseTileIterator
  use Grid_tile,           ONLY : Grid_tile_t
  use Grid_iterator,       ONLY : Grid_iterator_t
  use Timers_interface,    ONLY : Timers_start, Timers_stop
  use Driver_interface,    ONLY : Driver_getNstep
  use Stencils_interface,  ONLY : Stencils_advectWeno2d, Stencils_advectWeno3d
  use IncompNS_data

!------------------------------------------------------------------------------------------
  implicit none
  include "Flashx_mpi.h"
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
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
  !
  CALL SYSTEM_CLOCK(TA(1),count_rate)
  !
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
     !
#if NDIM == 3
     call tileDesc%getDataPtr(facezData, FACEZ)
     ! compute RHS of momentum equation
     call Stencils_advectWeno3d(facexData(RHDS_FACE_VAR,:,:,:), &
                                facexData(VELC_FACE_VAR,:,:,:), &
                                facexData(VELC_FACE_VAR,:,:,:), &
                                faceyData(VELC_FACE_VAR,:,:,:), &
                                facezData(VELC_FACE_VAR,:,:,:), &
                                del(DIR_X), del(DIR_Y), del(DIR_Z),&
                                GRID_ILO, GRID_IHI+1,&
                                GRID_JLO, GRID_JHI,& 
                                GRID_KLO, GRID_KHI,&
                                center=.false.,facex=.true.,facey=.false.,facez=.false.)

     call Stencils_advectWeno3d(faceyData(RHDS_FACE_VAR,:,:,:), &
                                faceyData(VELC_FACE_VAR,:,:,:), &
                                facexData(VELC_FACE_VAR,:,:,:), &
                                faceyData(VELC_FACE_VAR,:,:,:), &
                                facezData(VELC_FACE_VAR,:,:,:), &
                                del(DIR_X), del(DIR_Y), del(DIR_Z),&
                                GRID_ILO, GRID_IHI,&
                                GRID_JLO, GRID_JHI+1,& 
                                GRID_KLO, GRID_KHI,&
                                center=.false.,facex=.false.,facey=.true.,facez=.false.)

     call Stencils_advectWeno3d(facezData(RHDS_FACE_VAR,:,:,:), &
                                facezData(VELC_FACE_VAR,:,:,:), &
                                facexData(VELC_FACE_VAR,:,:,:), &
                                faceyData(VELC_FACE_VAR,:,:,:), &
                                facezData(VELC_FACE_VAR,:,:,:), &
                                del(DIR_X), del(DIR_Y), del(DIR_Z),&
                                GRID_ILO, GRID_IHI,&
                                GRID_JLO, GRID_JHI,& 
                                GRID_KLO, GRID_KHI+1,&
                                center=.false.,facex=.false.,facey=.false.,facez=.true.)
#elif NDIM ==2
     ! compute RHS of momentum equation
     call Stencils_advectWeno2d(facexData(RHDS_FACE_VAR,:,:,:), &
                                facexData(VELC_FACE_VAR,:,:,:), &
                                facexData(VELC_FACE_VAR,:,:,:), &
                                faceyData(VELC_FACE_VAR,:,:,:), &
                                del(DIR_X), &
                                del(DIR_Y), &
                                GRID_ILO, GRID_IHI+1,&
                                GRID_JLO, GRID_JHI,& 
                                center=.false.,facex=.true.,facey=.false.)

     call Stencils_advectWeno2d(faceyData(RHDS_FACE_VAR,:,:,:), &
                                faceyData(VELC_FACE_VAR,:,:,:), &
                                facexData(VELC_FACE_VAR,:,:,:), &
                                faceyData(VELC_FACE_VAR,:,:,:), &
                                del(DIR_X), &
                                del(DIR_Y), &
                                GRID_ILO, GRID_IHI,&
                                GRID_JLO, GRID_JHI+1,& 
                                center=.false.,facex=.false.,facey=.true.)

#endif
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

  CALL SYSTEM_CLOCK(TA(2),count_rate)
  ET=REAL(TA(2)-TA(1))/count_rate
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total IncompNS Advection Time =',ET

  return
end subroutine IncompNS_advection
