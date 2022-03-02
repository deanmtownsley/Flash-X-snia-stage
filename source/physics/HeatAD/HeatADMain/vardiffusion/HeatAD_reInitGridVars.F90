!!****if* source/physics/HeatAD/HeatADMain/vardiffusion/HeatAD_reInitGridVars
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

subroutine HeatAD_reInitGridVars()

  use Grid_interface,   ONLY : Grid_getTileIterator,Grid_releaseTileIterator
  use Grid_tile,        ONLY : Grid_tile_t
  use Grid_iterator,    ONLY : Grid_iterator_t
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_getNStep
  use HeatAD_data

!------------------------------------------------------------------------------------------
  implicit none
  include "Flashx_mpi.h"
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: solnData
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor

!------------------------------------------------------------------------------------------
  nullify(solnData)

  call Timers_start("HeatAD_reInitGridVars")

  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData,  CENTER)
     solnData(ALPH_VAR,:,:,:) = 1.
     solnData(RHST_VAR,:,:,:) = 0.
     solnData(TFRC_VAR,:,:,:) = 0.
     ! Release pointers:
     call tileDesc%releaseDataPtr(solnData,  CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)  

  call Timers_stop("HeatAD_reInitGridVars")

  return
end subroutine HeatAD_reInitGridVars
