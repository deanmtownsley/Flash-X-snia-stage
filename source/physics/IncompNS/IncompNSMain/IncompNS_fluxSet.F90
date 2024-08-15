!!****if* source/physics/IncompNS/IncompNSMain/constDens/IncompNS_fluxSet
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
!!REORDER(4): flux[xyz]Data

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

subroutine IncompNS_fluxSet(tileDesc)

   use Grid_interface, ONLY: Grid_putFluxData
   use Grid_tile, ONLY: Grid_tile_t
   use Timers_interface, ONLY: Timers_start, Timers_stop

   implicit none
   type(Grid_tile_t), INTENT(IN) :: tileDesc

   real :: del(MDIM)
   integer :: lo(3), hi(3)
#if NDIM < MDIM
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, fluxxData, fluxyData, fluxzData
#else
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData, fluxxData, fluxyData, fluxzData
#endif
!----------------------------------------------------------------------------------------
#if NDIM < MDIM
   nullify (solnData, facexData, faceyData, fluxxData, fluxyData, fluxzData)
#else
   nullify (solnData, facexData, faceyData, facezData, fluxxData, fluxyData, fluxzData)
#endif

   call Timers_start("IncompNS_fluxSet")

   call tileDesc%deltas(del)
#if NDIM == 2
   del(DIR_Z) = 1
#endif

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%getDataPtr(facexData, FACEX)
   call tileDesc%getDataPtr(faceyData, FACEY)
   call tileDesc%getDataPtr(fluxxData, FLUXX)
   call tileDesc%getDataPtr(fluxyDAta, FLUXY)
#if NDIM == 3
   call tileDesc%getDataPtr(facezData, FACEZ)
#endif
   call tileDesc%getDataPtr(fluxzData, FLUXZ)

   lo(1:MDIM) = tileDesc%limits(LOW, 1:MDIM)
   hi(1:MDIM) = tileDesc%limits(HIGH, 1:MDIM)

#ifdef FLASH_GRID_AMREX
   fluxxData(MOMT_FLUX, :, :, :) = facexData(VELC_FACE_VAR, lo(1):hi(1)+1, lo(2):hi(2), lo(3):hi(3))
   fluxyData(MOMT_FLUX, :, :, :) = faceyData(VELC_FACE_VAR, lo(1):hi(1), lo(2):hi(2)+1, lo(3):hi(3))
#if NDIM==3
   fluxzData(MOMT_FLUX, :, :, :) = facezData(VELC_FACE_VAR, lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)+1)
#endif

   call Grid_putFluxData(tileDesc, fluxxData, fluxyData, fluxzData, lo)
#endif

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)
   call tileDesc%releaseDataPtr(fluxxData, FLUXX)
   call tileDesc%releaseDataPtr(fluxyData, FLUXY)
#if NDIM ==3
   call tileDesc%releaseDataPtr(facezData, FACEZ)
   call tileDesc%releaseDataPtr(fluxzData, FLUXZ)
#endif

   call Timers_stop("IncompNS_fluxSet")

end subroutine IncompNS_fluxSet
