!!****if* source/physics/IncompNS/IncompNSMain/varDens/IncompNS_correctorUpdate
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
!!REORDER(4): fluxBuf[XYZ]

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

subroutine IncompNS_correctorUpdate(tileDesc, fluxBufX, fluxBufY, fluxBufZ, lo)

   use Grid_tile, ONLY: Grid_tile_t
   use Timers_interface, ONLY: Timers_start, Timers_stop

   implicit none
   type(Grid_tile_t), INTENT(IN) :: tileDesc
   integer, intent(in) :: lo(3)
   real, intent(in), dimension(1:, lo(1):, lo(2):, lo(3):) :: fluxBufX, fluxBufY, fluxBufZ

   real :: del(MDIM)
   integer :: sx, sy, sz, ex, ey, ez
   integer :: level
#if NDIM < MDIM
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData
   real, dimension(NFACE_VARS, 1, 1, 1) :: facezData
#else
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
#endif
!----------------------------------------------------------------------------------------
#if NDIM < MDIM
   nullify (solnData, facexData, faceyData)
#else
   nullify (solnData, facexData, faceyData, facezData)
#endif

   call Timers_start("IncompNS_correctorUpdate")

   call tileDesc%deltas(del)

#if NDIM == 2
   del(DIR_Z) = 1
#endif

   ! Positions in face arrays where flux vars have been stored
   sx = tileDesc%limits(LOW, IAXIS)
   sy = tileDesc%limits(LOW, JAXIS)
   sz = 1

   ex = tileDesc%limits(HIGH, IAXIS)
   ey = tileDesc%limits(HIGH, JAXIS)
   ez = 1

#if NDIM == 3
   sz = tileDesc%limits(LOW, KAXIS)
   ez = tileDesc%limits(HIGH, KAXIS)
#endif

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%getDataPtr(facexData, FACEX)
   call tileDesc%getDataPtr(faceyData, FACEY)

#if NDIM == 3
   call tileDesc%getDataPtr(facezData, FACEZ)
#endif

   ! update interface velocities with flux variables
   ! X direction:
   facexData(VELC_FACE_VAR, sx, sy:ey, sz:ez) = &
      fluxBufX(VELC_FLUX, sx, sy:ey, sz:ez)/(del(DIR_Y)*del(DIR_Z))  ! dy*dz
   facexData(VELC_FACE_VAR, ex + 1, sy:ey, sz:ez) = &
      fluxBufX(VELC_FLUX, ex + 1, sy:ey, sz:ez)/(del(DIR_Y)*del(DIR_Z))

   ! Y direction:
   faceyData(VELC_FACE_VAR, sx:ex, sy, sz:ez) = &
      fluxBufY(VELC_FLUX, sx:ex, sy, sz:ez)/(del(DIR_X)*del(DIR_Z))  ! dx*dz
   faceyData(VELC_FACE_VAR, sx:ex, ey + 1, sz:ez) = &
      fluxBufY(VELC_FLUX, sx:ex, ey + 1, sz:ez)/(del(DIR_X)*del(DIR_Z))

#if NDIM == 3
   ! Z direction:
   facezData(VELC_FACE_VAR, sx:ex, sy:ey, sz) = &
      fluxBufZ(VELC_FLUX, sx:ex, sy:ey, sz)/(del(DIR_X)*del(DIR_Y))  ! dx*dy
   facezData(VELC_FACE_VAR, sx:ex, sy:ey, ez + 1) = &
      fluxBufZ(VELC_FLUX, sx:ex, sy:ey, ez + 1)/(del(DIR_X)*del(DIR_Y))
#endif

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)

#if NDIM ==3
   call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

   call Timers_stop("IncompNS_correctorUpdate")

end subroutine IncompNS_correctorUpdate
