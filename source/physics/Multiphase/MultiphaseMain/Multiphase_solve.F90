!!***if* source/physics/Multiphase/MultiphaseMain/Multiphase_solve
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

subroutine Multiphase_solve(dt)

   use Multiphase_data
   use Timers_interface, ONLY: Timers_start, Timers_stop
   use Driver_interface, ONLY: Driver_getNStep
   use Grid_interface, ONLY: Grid_getTileIterator, Grid_releaseTileIterator
   use Grid_tile, ONLY: Grid_tile_t
   use Grid_iterator, ONLY: Grid_iterator_t
   use Stencils_interface, ONLY: Stencils_integrateEuler

   implicit none
   include "Flashx_mpi.h"
   !----------Arugments List---------------
   real, INTENT(IN) :: dt

!-----------------------------------------------------------------------------------------
   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:, :, :, :) :: solnData
   real del(MDIM)
   type(Grid_tile_t) :: tileDesc
   type(Grid_iterator_t) :: itor

!-----------------------------------------------------------------------------------------
   nullify (solnData)

   call Timers_start("Multiphase_solve")

   call Grid_getTileIterator(itor, nodetype=LEAF)
   do while (itor%isValid())
      call itor%currentTile(tileDesc)
      call tileDesc%getDataPtr(solnData, CENTER)

      call Stencils_integrateEuler(solnData(DFUN_VAR, :, :, :), &
                                   solnData(RDFN_VAR, :, :, :), &
                                   dt, &
                                   GRID_ILO, GRID_IHI, &
                                   GRID_JLO, GRID_JHI, &
                                   GRID_KLO, GRID_KHI, &
                                   iSource=solnData(DFRC_VAR, :, :, :))

      ! Release pointers:
      call tileDesc%releaseDataPtr(solnData, CENTER)
      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)

   call Timers_stop("Multiphase_solve")

   return
end subroutine Multiphase_solve
