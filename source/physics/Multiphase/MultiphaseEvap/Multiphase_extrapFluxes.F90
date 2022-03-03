!!***if* source/physics/Multiphase/MultiphaseEvap/Multiphase_extrapFluxes
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
#include "Multiphase.h"
#include "Simulation.h"

subroutine Multiphase_extrapFluxes(tileDesc, iteration)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t
   use Stencils_interface, ONLY: Stencils_cnt_advectUpwind2d, Stencils_cnt_advectUpwind3d
   use mph_evapInterface, ONLY: mph_phasedFluxes

!------------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   integer, intent(in) :: iteration
   type(Grid_tile_t), intent(in) :: tileDesc

   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:, :, :, :) :: solnData
   real del(MDIM)
!------------------------------------------------------------------------------------------------
   nullify (solnData)

   call Timers_start("Multiphase_extrapFluxes")

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%deltas(del)

#if NDIM < MDIM

   solnData(MFLX_VAR, :, :, :) = 0.0

   call Stencils_cnt_advectUpwind2d(solnData(MFLX_VAR, :, :, :), &
                                    solnData(HFLQ_VAR, :, :, :), &
                                    solnData(NRMX_VAR, :, :, :), &
                                    solnData(NRMY_VAR, :, :, :), &
                                    del(IAXIS), del(JAXIS), &
                                    GRID_ILO, GRID_IHI, &
                                    GRID_JLO, GRID_JHI)

   call mph_phasedFluxes(solnData(HFLQ_VAR, :, :, :), &
                         solnData(MFLX_VAR, :, :, :), &
                         solnData(DFUN_VAR, :, :, :), &
                         0.5*del(IAXIS), &
                         GRID_ILO, GRID_IHI, &
                         GRID_JLO, GRID_JHI, &
                         GRID_KLO, GRID_KHI)

   solnData(MFLX_VAR, :, :, :) = 0.0

   call Stencils_cnt_advectUpwind2d(solnData(MFLX_VAR, :, :, :), &
                                    solnData(HFGS_VAR, :, :, :), &
                                    -solnData(NRMX_VAR, :, :, :), &
                                    -solnData(NRMY_VAR, :, :, :), &
                                    del(IAXIS), del(JAXIS), &
                                    GRID_ILO, GRID_IHI, &
                                    GRID_JLO, GRID_JHI)

   call mph_phasedFluxes(solnData(HFGS_VAR, :, :, :), &
                         solnData(MFLX_VAR, :, :, :), &
                         -solnData(DFUN_VAR, :, :, :), &
                         0.5*del(IAXIS), &
                         GRID_ILO, GRID_IHI, &
                         GRID_JLO, GRID_JHI, &
                         GRID_KLO, GRID_KHI)

#else

   solnData(MFLX_VAR, :, :, :) = 0.0

   call Stencils_cnt_advectUpwind3d(solnData(MFLX_VAR, :, :, :), &
                                    solnData(HFLQ_VAR, :, :, :), &
                                    solnData(NRMX_VAR, :, :, :), &
                                    solnData(NRMY_VAR, :, :, :), &
                                    solnData(NRMZ_VAR, :, :, :), &
                                    del(IAXIS), del(JAXIS), del(KAXIS), &
                                    GRID_ILO, GRID_IHI, &
                                    GRID_JLO, GRID_JHI, &
                                    GRID_KLO, GRID_KHI)

   call mph_phasedFluxes(solnData(HFLQ_VAR, :, :, :), &
                         solnData(MFLX_VAR, :, :, :), &
                         solnData(DFUN_VAR, :, :, :), &
                         0.5*del(IAXIS), &
                         GRID_ILO, GRID_IHI, &
                         GRID_JLO, GRID_JHI, &
                         GRID_KLO, GRID_KHI)

   solnData(MFLX_VAR, :, :, :) = 0.0

   call Stencils_cnt_advectUpwind3d(solnData(MFLX_VAR, :, :, :), &
                                    solnData(HFGS_VAR, :, :, :), &
                                    -solnData(NRMX_VAR, :, :, :), &
                                    -solnData(NRMY_VAR, :, :, :), &
                                    -solnData(NRMZ_VAR, :, :, :), &
                                    del(IAXIS), del(JAXIS), del(KAXIS), &
                                    GRID_ILO, GRID_IHI, &
                                    GRID_JLO, GRID_JHI, &
                                    GRID_KLO, GRID_KHI)

   call mph_phasedFluxes(solnData(HFGS_VAR, :, :, :), &
                         solnData(MFLX_VAR, :, :, :), &
                         -solnData(DFUN_VAR, :, :, :), &
                         0.5*del(IAXIS), &
                         GRID_ILO, GRID_IHI, &
                         GRID_JLO, GRID_JHI, &
                         GRID_KLO, GRID_KHI)

#endif

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)

   call Timers_stop("Multiphase_extrapFluxes")

end subroutine Multiphase_extrapFluxes
