!!***if* source/physics/Multiphase/MultiphaseMain/Multiphase_setPressureJumps
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
!!REORDER(4): solnData, face[xyz]Data

#include "constants.h"
#include "Multiphase.h"
#include "Simulation.h"

subroutine Multiphase_setPressureJumps(tileDesc)

   use Multiphase_data
   use mph_interface, ONLY: mph_setWeberJumps2d, mph_setWeberJumps3d
   use mph_evapInterface, ONLY: mph_setEvapJumps2d, mph_setEvapJumps3d
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t
   use Stencils_interface, ONLY: Stencils_lsNormals2d, Stencils_lsNormals3d

!----------------------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   integer :: ierr, i, j, k
   real del(MDIM)
!---------------------------------------------------------------------------------------------------------
   nullify (solnData, facexData, faceyData, facezData)

   call Timers_start("Multiphase_setPressureJumps")

   !-------------------------------------------------
   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%getDataPtr(facexData, FACEX)
   call tileDesc%getDataPtr(faceyData, FACEY)
#if NDIM == MDIM
   call tileDesc%getDataPtr(facezData, FACEZ)
#endif
   call tileDesc%deltas(del)

#if NDIM < MDIM
   call mph_setWeberJumps2d(solnData(DFUN_VAR, :, :, :), &
                            solnData(CURV_VAR, :, :, :), &
                            solnData(PFUN_VAR, :, :, :), &
                            facexData(mph_iJumpVar, :, :, :), &
                            faceyData(mph_iJumpVar, :, :, :), &
                            del(DIR_X), del(DIR_Y), &
                            mph_invWeber, mph_rhoGas, &
                            GRID_ILO, GRID_IHI, &
                            GRID_JLO, GRID_JHI, tol=mph_presTol)

#ifdef MULTIPHASE_EVAPORATION
   call mph_setEvapJumps2d(solnData(DFUN_VAR, :, :, :), &
                           solnData(PFUN_VAR, :, :, :), &
                           facexData(mph_iJumpVar, :, :, :), &
                           faceyData(mph_iJumpVar, :, :, :), &
                           solnData(MFLX_VAR, :, :, :), mph_rhoGas, &
                           del(DIR_X), del(DIR_Y), &
                           GRID_ILO, GRID_IHI, &
                           GRID_JLO, GRID_JHI, tol=mph_presTol)
#endif

#else

   call mph_setWeberJumps3d(solnData(DFUN_VAR, :, :, :), &
                            solnData(CURV_VAR, :, :, :), &
                            solnData(PFUN_VAR, :, :, :), &
                            facexData(mph_iJumpVar, :, :, :), &
                            faceyData(mph_iJumpVar, :, :, :), &
                            facezData(mph_iJumpVar, :, :, :), &
                            del(DIR_X), del(DIR_Y), del(DIR_Z), &
                            mph_invWeber, mph_rhoGas, &
                            GRID_ILO, GRID_IHI, &
                            GRID_JLO, GRID_JHI, &
                            GRID_KLO, GRID_KHI, tol=mph_presTol)

#ifdef MULTIPHASE_EVAPORATION
   call mph_setEvapJumps3d(solnData(DFUN_VAR, :, :, :), &
                           solnData(PFUN_VAR, :, :, :), &
                           facexData(mph_iJumpVar, :, :, :), &
                           faceyData(mph_iJumpVar, :, :, :), &
                           facezData(mph_iJumpVar, :, :, :), &
                           solnData(MFLX_VAR, :, :, :), mph_rhoGas, &
                           del(DIR_X), del(DIR_Y), del(DIR_Z), &
                           GRID_ILO, GRID_IHI, &
                           GRID_JLO, GRID_JHI, &
                           GRID_KLO, GRID_KHI, tol=mph_presTol)
#endif

#endif
   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)
   call tileDesc%releaseDataPtr(facexData, FACEX)
   call tileDesc%releaseDataPtr(faceyData, FACEY)
#if NDIM == MDIM
   call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

   call Timers_stop("Multiphase_setPressureJumps")

   return

end subroutine Multiphase_setPressureJumps
