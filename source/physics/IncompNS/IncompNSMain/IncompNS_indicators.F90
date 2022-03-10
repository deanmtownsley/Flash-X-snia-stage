!!****if* source/physics/IncompNS/IncompNSMain/IncompNS_indicators
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

subroutine IncompNS_indicators()

   use Grid_interface, ONLY: Grid_getTileIterator, Grid_releaseTileIterator
   use Grid_tile, ONLY: Grid_tile_t
   use Grid_iterator, ONLY: Grid_iterator_t
   use ins_interface, ONLY: ins_indicators
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use IncompNS_data

!------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   integer, dimension(2, MDIM) :: stnLimits = 1
   real :: vecminaux(5), vecmaxaux(5), vecmin(5), vecmax(5)
#if NDIM < MDIM
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData
   real, dimension(NFACE_VARS, 1, 1, 1) :: facezData
#else
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
#endif
   integer :: NStep, ierr
   type(Grid_tile_t) :: tileDesc
   type(Grid_iterator_t) :: itor

!------------------------------------------------------------------------------------------
#if NDIM < MDIM
   nullify (solnData, facexData, faceyData)
#else
   nullify (solnData, facexData, faceyData, facezData)
#endif

   vecmaxaux = -10.**(10.)
   vecminaux = 10.**(10.)
   call Grid_getTileIterator(itor, nodetype=LEAF)
   do while (itor%isValid())
      call itor%currentTile(tileDesc)

      stnLimits(LOW, 1:NDIM) = tileDesc%limits(LOW, 1:NDIM) - tileDesc%blkLimitsGC(LOW, 1:NDIM) + 1
      stnLimits(HIGH, 1:NDIM) = tileDesc%limits(HIGH, 1:NDIM) - tileDesc%blkLimitsGC(LOW, 1:NDIM) + 1

      call tileDesc%getDataPtr(solnData, CENTER)
      call tileDesc%getDataPtr(facexData, FACEX)
      call tileDesc%getDataPtr(faceyData, FACEY)
#if NDIM == 3
      call tileDesc%getDataPtr(facezData, FACEZ)
#endif
      ! compute indicators
      call ins_indicators(facexData(VELC_FACE_VAR, :, :, :), &
                          faceyData(VELC_FACE_VAR, :, :, :), &
                          facezData(VELC_FACE_VAR, :, :, :), &
                          solnData(PRES_VAR, :, :, :), &
                          solnData(DUST_VAR, :, :, :), &
                          stnLimits(LOW, IAXIS), stnLimits(HIGH, IAXIS), &
                          stnLimits(LOW, JAXIS), stnLimits(HIGH, JAXIS), &
                          stnLimits(LOW, KAXIS), stnLimits(HIGH, KAXIS), &
                          vecminaux, vecmaxaux)

      ! Release pointers:
      call tileDesc%releaseDataPtr(solnData, CENTER)
      call tileDesc%releaseDataPtr(facexData, FACEX)
      call tileDesc%releaseDataPtr(faceyData, FACEY)
#if NDIM ==3
      call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif
      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)

   call MPI_Allreduce(vecmaxaux, vecmax, 5, FLASH_REAL, &
                      MPI_MAX, ins_meshComm, ierr)

   call MPI_Allreduce(vecminaux, vecmin, 5, FLASH_REAL, &
                      MPI_MIN, ins_meshComm, ierr)

   if (ins_meshMe .eq. MASTER_PE) then
      write (*, *) ' '
      write (*, '(A24,2g14.6)') ' Min , Max  U =', vecmin(2), vecmax(2) !minu,maxu
      write (*, '(A24,2g14.6)') ' Min , Max  V =', vecmin(3), vecmax(3) !minv,maxv
#if NDIM == 3
      write (*, '(A24,2g14.6)') ' Min , Max  W =', vecmin(4), vecmax(4) !minw,maxw
#endif
      write (*, '(A24,2g14.6)') ' Min , Max  P =', vecmin(5), vecmax(5) !minp,maxp
      write (*, '(A24,2g14.6)') ' Min , Max  Divergence =', vecmin(1), vecmax(1) !mndivv,mxdivv
   end if
   !----------------------------------------------------------------------------------------------------

   ins_mindiv = vecmin(1)
   ins_maxdiv = vecmax(1)

   return
end subroutine IncompNS_indicators
