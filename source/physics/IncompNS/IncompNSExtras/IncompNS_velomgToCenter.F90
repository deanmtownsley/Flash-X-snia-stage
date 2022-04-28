!!****if* source/physics/IncompNS/IncompNSExtras/IncompNS_velomgToCenter
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
!!
!!  IncompNS_velomg2center
!!
!! SYNOPSIS
!!
!!  call IncompNS_velomg2center(integer(in) :: blockList,
!!                              integer(in) :: blockCount)
!!
!! DESCRIPTION
!!
!!   Compute cell-centered (or volume-averaged) versions of some
!!   quantities, by interpolation or averaging from face-centered
!!   versions.
!!
!! ARGUMENTS
!!
!!   blockList : list of block IDs
!!
!!   blockCount : number of block IDs in the list
!!
!!
!! NOTES
!!
!!  To write cell centered velocities and vorticity velx,vely,velz,omgx,omgy,omgz,
!!  add REQUIRES physics/IncompNS/IncompNSExtras to Config.
!!
!!***
!!REORDER(4): face[xyz]Data
!!REORDER(4): solnData

subroutine IncompNS_velomgToCenter()

  use Grid_interface,   ONLY : Grid_fillGuardCells, Grid_getTileIterator, &
                               Grid_releaseTileIterator, Grid_getNumBlksFromType, &
                               Grid_getCellVolumes
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t
  use ins_extrasInterface, ONLY: ins_velToCenter, ins_omgToCenter

  implicit none
#include "constants.h"
#include "Simulation.h"

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  real :: del(MDIM)
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)

  nullify(solnData,facexData,faceyData,facezData)

#if NDIM < MDIM
allocate(facezData(VELC_FACE_VAR,1,1,1))
#endif

  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())

     call itor%currentTile(tileDesc)

     call tileDesc%deltas(del)

     blkLimits   = tileDesc%limits
     blkLimitsGC = tileDesc%blkLimitsGC

     call tileDesc%getDataPtr(solnData,  CENTER)
     call tileDesc%getDataPtr(facexData, FACEX)
     call tileDesc%getDataPtr(faceyData, FACEY)

#if NDIM == 3
     call tileDesc%getDataPtr(facezData, FACEZ)
#endif

     call ins_velToCenter(facexData(VELC_FACE_VAR,:,:,:),&
                          faceyData(VELC_FACE_VAR,:,:,:),&
                          facezData(VELC_FACE_VAR,:,:,:),&
                          solnData(VELX_VAR,:,:,:),&
                          solnData(VELY_VAR,:,:,:),&
                          solnData(VELZ_VAR,:,:,:),&
                          GRID_ILO_GC,GRID_IHI_GC,&
                          GRID_JLO_GC,GRID_JHI_GC,&
                          GRID_KLO_GC,GRID_KHI_GC)

     call ins_omgToCenter(facexData(VELC_FACE_VAR,:,:,:),&
                          faceyData(VELC_FACE_VAR,:,:,:),&
                          facezData(VELC_FACE_VAR,:,:,:),&
                          solnData(OMGM_VAR,:,:,:),&
                          GRID_ILO_GC,GRID_IHI_GC,&
                          GRID_JLO_GC,GRID_JHI_GC,&
                          GRID_KLO_GC,GRID_KHI_GC,&
                          del(IAXIS),del(JAXIS),del(KAXIS))

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

#if NDIM < MDIM
deallocate(facezData)
#endif

  return

end subroutine IncompNS_velomgToCenter
