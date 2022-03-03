!!***if* source/physics/Multiphase/MultiphaseEvap/Multiphase_setMassFlux
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

subroutine Multiphase_setMassFlux(tileDesc)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_tile, ONLY: Grid_tile_t
   use Stencils_interface, ONLY: Stencils_cnt_advectUpwind2d, Stencils_cnt_advectUpwind3d
   use mph_evapInterface, ONLY: mph_phasedFluxes

!------------------------------------------------------------------------------------------------
   implicit none
   include "Flashx_mpi.h"
   type(Grid_tile_t), intent(in) :: tileDesc

   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:, :, :, :) :: solnData
   real del(MDIM)

!------------------------------------------------------------------------------------------------
   nullify (solnData)

   call Timers_start("Multiphase_setMassFlux")

   call tileDesc%getDataPtr(solnData, CENTER)
   call tileDesc%deltas(del)

   solnData(MFLX_VAR, :, :, :) = (mph_Stefan*mph_invReynolds/mph_Prandtl)* &
                                 (solnData(HFLQ_VAR, :, :, :) + mph_thcoGas*solnData(HFGS_VAR, :, :, :))

   ! Release pointers:
   call tileDesc%releaseDataPtr(solnData, CENTER)

   call Timers_stop("Multiphase_setMassFlux")

end subroutine Multiphase_setMassFlux
