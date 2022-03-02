!!****if* source/physics/IncompNS/IncompNSMain/vardens/IncompNS_corrector
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
!!
!!***
!!REORDER(4): face[xyz]Data
!!REORDER(4): solnData

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

subroutine IncompNS_corrector(dt)

  use Grid_interface,   ONLY : Grid_getTileIterator, &
                               Grid_releaseTileIterator
  use Grid_tile,        ONLY : Grid_tile_t
  use Grid_iterator,    ONLY : Grid_iterator_t
  use ins_interface,    ONLY : ins_corrector_constdens
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_getNStep
  use IncompNS_data

  implicit none
  include "Flashx_mpi.h"
  !-----Argument List-----
  real,    INTENT(IN) :: dt

!------------------------------------------------------------------------------------------
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
#if NDIM < MDIM
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData
  real, dimension(NFACE_VARS,1,1,1) :: facezData
#else
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
#endif
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

  call Timers_start("IncompNS_corrector")

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
     call ins_corrector_constdens(facexData(VELC_FACE_VAR,:,:,:),&
                                  faceyData(VELC_FACE_VAR,:,:,:),&
                                  facezData(VELC_FACE_VAR,:,:,:),&
                                  facexData(PGN1_FACE_VAR,:,:,:),&
                                  faceyData(PGN1_FACE_VAR,:,:,:),&
                                  facezData(PGN1_FACE_VAR,:,:,:),&                           
                                  solnData(PRES_VAR,:,:,:),&
                                  GRID_ILO,GRID_IHI,&
                                  GRID_JLO,GRID_JHI,&
                                  GRID_KLO,GRID_KHI,&
                                  dt,del(DIR_X),del(DIR_Y),del(DIR_Z))

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

  call Timers_stop("IncompNS_corrector")

  return
end subroutine IncompNS_corrector
