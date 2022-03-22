!!****if* source/physics/IncompNS/IncompNSMain/varDens/IncompNS_advection
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

subroutine IncompNS_advection(tileDesc)

   use Grid_tile, ONLY: Grid_tile_t
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNstep
   use Stencils_interface, ONLY: Stencils_advectWeno2d, Stencils_advectWeno3d
   use IncompNS_data

!------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
#if NDIM < MDIM
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData
   real, dimension(NFACE_VARS, 1, 1, 1) :: facezData
#else
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
#endif
   real del(MDIM)
   integer :: NStep

!------------------------------------------------------------------------------------------
#if NDIM < MDIM
   nullify (solnData, facexData, faceyData)
#else
   nullify (solnData, facexData, faceyData, facezData)
#endif
   !
   call Timers_start("IncompNS_advection")
   !
   blkLimits = tileDesc%limits
   blkLimitsGC = tileDesc%blkLimitsGC

   call tileDesc%deltas(del)
   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%getDataPtr(facexData, FACEX)
   call tileDesc%getDataPtr(faceyData, FACEY)
   !
#if NDIM == 3
   call tileDesc%getDataPtr(facezData, FACEZ)
   ! compute RHS of momentum equation
   call Stencils_advectWeno3d(facexData(HVN0_FACE_VAR, :, :, :), &
                              facexData(VELC_FACE_VAR, :, :, :), &
                              facexData(VELC_FACE_VAR, :, :, :), &
                              faceyData(VELC_FACE_VAR, :, :, :), &
                              facezData(VELC_FACE_VAR, :, :, :), &
                              del(DIR_X), del(DIR_Y), del(DIR_Z), &
                              GRID_ILO, GRID_IHI + 1, &
                              GRID_JLO, GRID_JHI, &
                              GRID_KLO, GRID_KHI, &
                              center=0, facex=1, facey=0, facez=0)

   call Stencils_advectWeno3d(faceyData(HVN0_FACE_VAR, :, :, :), &
                              faceyData(VELC_FACE_VAR, :, :, :), &
                              facexData(VELC_FACE_VAR, :, :, :), &
                              faceyData(VELC_FACE_VAR, :, :, :), &
                              facezData(VELC_FACE_VAR, :, :, :), &
                              del(DIR_X), del(DIR_Y), del(DIR_Z), &
                              GRID_ILO, GRID_IHI, &
                              GRID_JLO, GRID_JHI + 1, &
                              GRID_KLO, GRID_KHI, &
                              center=0, facex=0, facey=1, facez=0)

   call Stencils_advectWeno3d(facezData(HVN0_FACE_VAR, :, :, :), &
                              facezData(VELC_FACE_VAR, :, :, :), &
                              facexData(VELC_FACE_VAR, :, :, :), &
                              faceyData(VELC_FACE_VAR, :, :, :), &
                              facezData(VELC_FACE_VAR, :, :, :), &
                              del(DIR_X), del(DIR_Y), del(DIR_Z), &
                              GRID_ILO, GRID_IHI, &
                              GRID_JLO, GRID_JHI, &
                              GRID_KLO, GRID_KHI + 1, &
                              center=0, facex=0, facey=0, facez=1)
#elif NDIM ==2
   ! compute RHS of momentum equation
   call Stencils_advectWeno2d(facexData(HVN0_FACE_VAR, :, :, :), &
                              facexData(VELC_FACE_VAR, :, :, :), &
                              facexData(VELC_FACE_VAR, :, :, :), &
                              faceyData(VELC_FACE_VAR, :, :, :), &
                              del(DIR_X), &
                              del(DIR_Y), &
                              GRID_ILO, GRID_IHI + 1, &
                              GRID_JLO, GRID_JHI, &
                              center=0, facex=1, facey=0)

   call Stencils_advectWeno2d(faceyData(HVN0_FACE_VAR, :, :, :), &
                              faceyData(VELC_FACE_VAR, :, :, :), &
                              facexData(VELC_FACE_VAR, :, :, :), &
                              faceyData(VELC_FACE_VAR, :, :, :), &
                              del(DIR_X), &
                              del(DIR_Y), &
                              GRID_ILO, GRID_IHI, &
                              GRID_JLO, GRID_JHI + 1, &
                              center=0, facex=0, facey=1)

#endif
   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)

#if NDIM ==3
   call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

   call Timers_stop("IncompNS_advection")

   return
end subroutine IncompNS_advection
